/***********************************************************************
 * @file		ParallelGenomicChunks.h
 * @author		Zilong Li
 * @brief		get chunks in parallel for pca stream
 * Copyright (C) 2025. The use of this code is governed by the LICENSE file.
 ***********************************************************************/

#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
#include <mutex>
#include <condition_variable>
#include <thread>
#include <atomic>
#include <future>
#include <functional>
#include <memory>
#include <chrono>
#include <type_traits>

#include <RcppEigen.h>

#include "GenomicDataStream.h"
#include "GenomicDataStream_virtual.h"
#include "threadpool.hpp"

// Thread-safe chunk cache for coordinating between reader and worker threads
template<typename T>
class ChunkCache {
public:
  ChunkCache(size_t maxQueueSize) : maxSize(maxQueueSize), done(false) {}

  void addChunk(std::shared_ptr<gds::DataChunk<T>> chunk, size_t chunkIndex) {
    std::unique_lock<std::mutex> lock(mutex);
    // Wait until there's room in the queue
    notFull.wait(lock, [this] { return queue.size() < maxSize || done; });

    if (done)
      return;

    queue.push(std::make_pair(chunkIndex, chunk));
    lock.unlock();
    notEmpty.notify_one();
  }

  bool getChunk(std::shared_ptr<gds::DataChunk<T>> &chunk, size_t &chunkIndex) {
    std::unique_lock<std::mutex> lock(mutex);
    // Wait until there are items in the queue
    notEmpty.wait(lock, [this] { return !queue.empty() || done; });

    if (queue.empty() && done)
      return false;

    auto item = queue.front();
    queue.pop();
    chunkIndex = item.first;
    chunk = item.second;

    lock.unlock();
    notFull.notify_one();
    return true;
  }

  void setDone() {
    std::lock_guard<std::mutex> lock(mutex);
    done = true;
    notEmpty.notify_all();
    notFull.notify_all();
  }

private:
  std::queue<std::pair<size_t, std::shared_ptr<gds::DataChunk<T>>>> queue;
  std::mutex mutex;
  std::condition_variable notEmpty;
  std::condition_variable notFull;
  size_t maxSize;
  bool done;
};

// Manages the overall file processing with a dedicated reader thread
template<typename T>
class GenomicStreamProcessor {
public:
  GenomicStreamProcessor(const gds::Param &param, size_t numThreads)
      : chunkSize(param.chunkSize), numChunks(param.regions.size()),
        pool(numThreads),
        cache(std::min(numThreads * 2,
                       size_t(16))), // Cache size: 2 chunks per thread, max 16
        stopReading(false) {

    reader = gds::createFileView(param);

    std::cout << "Processing in " << numChunks << " chunks\n";

    // Pre-allocate chunk objects for the buffer pool
    const size_t maxPoolSize = 16;
    for (size_t i = 0; i < maxPoolSize; ++i) {
      chunkPool.push(std::make_shared<gds::DataChunk<T>>()); // nsnps x nsamples
    }

  }

  void
  processFile(std::function<void(const gds::DataChunk<T> &, size_t)> processChunk) {
    // Start reader thread
    std::thread readerThread([this] { readerFunction(); });

    std::vector<std::future<void>> futures;
    size_t processedChunks = 0;

    // Process chunks as they become available
    while (processedChunks < numChunks) {
      std::shared_ptr<gds::DataChunk<T>> chunk;
      size_t chunkIndex;

      if (!cache.getChunk(chunk, chunkIndex)) {
        break; // No more chunks available
      }

      // Submit task to thread pool

      futures.emplace_back(pool.enqueue([this, chunk, chunkIndex, processChunk]{
        // Process the chunk
        processChunk(*chunk, chunkIndex);
        // Return chunk to pool
        {
          std::lock_guard<std::mutex> lock(poolMutex);
          chunkPool.push(chunk);
        }
        poolCondition.notify_one();
      }));
      processedChunks++;
    }

    // Wait for all processing to complete
    for (auto &future : futures) {
      future.get();
    }

    // Stop the reader thread
    {
      std::lock_guard<std::mutex> lock(readerMutex);
      stopReading = true;
    }
    poolCondition.notify_one();

    if (readerThread.joinable()) {
      readerThread.join();
    }
  }

private:
  void readerFunction() {
    for (size_t i = 0; i < numChunks; ++i) {
      // Check if we should stop
      {
        std::lock_guard<std::mutex> lock(readerMutex);
        if (stopReading)
          break;
      }

      // Get a chunk from the pool or create new one if none available
      std::shared_ptr<gds::DataChunk<T>> chunk;
      {
        std::unique_lock<std::mutex> lock(poolMutex);
        poolCondition.wait(
            lock, [this] { return !chunkPool.empty() || stopReading; });

        if (stopReading)
          break;

        chunk = chunkPool.front();
        chunkPool.pop();
      }

      // Read and process chunk
      int bytesRead = reader->getNextChunk(*chunk);

      if (bytesRead > 0) {
        // Add to processing cache
        cache.addChunk(chunk, i);
      } else {
        // Return chunk to pool if no data was read
        std::lock_guard<std::mutex> lock(poolMutex);
        chunkPool.push(chunk);
      }
    }

    // Signal that there are no more chunks
    cache.setDone();
  }

  std::unique_ptr<gds::GenomicDataStream> reader = nullptr;

  const size_t chunkSize;
  const size_t numChunks;

  ThreadPool pool;
  ChunkCache<T> cache;

  // Chunk pool for recycling buffers
  std::queue<std::shared_ptr<gds::DataChunk<T>>> chunkPool;
  std::mutex poolMutex;
  std::condition_variable poolCondition;

  // Reader thread control
  std::mutex readerMutex;
  bool stopReading;
};

int get_chunks_parallel_vcf(const std::string &file, const std::string &region) {

  std::string field = "GT";  // read dosage field
  std::string samples = "-"; // no samples filter
  double minVariance = NAN;    // retain features with var > minVariance
  int chunkSize = 10;        // each chunk will read 10 variants
  gds::Param param(file, region, samples, minVariance, chunkSize);
  param.setField(field);

  // Configuration
  size_t NUM_THREADS = std::thread::hardware_concurrency();
  std::cout << "number of available threads" << NUM_THREADS << std::endl;
  NUM_THREADS = NUM_THREADS > 6 ? 6 : NUM_THREADS;
  std::mutex printMutex;
  std::atomic<size_t> totalProcessed{0};

  try {
    GenomicStreamProcessor<Eigen::MatrixXd> processor(param, NUM_THREADS);

    processor.processFile([&](const gds::DataChunk<Eigen::MatrixXd> &chunk, size_t chunkIndex) {
      // Example processing logic
      auto X = chunk.getData();
      size_t size = X.size();

      // Calculate sum of X for test
      double chunkSum = X.sum();

      // Thread-safe update
      {
        std::lock_guard<std::mutex> lock(printMutex);
        totalProcessed += size;
        std::cout << "Chunk " << chunkIndex << " processed: " << size
                  << " bytes, sum: " << chunkSum << std::endl;
      }
    });

  } catch (const std::exception &e) {
    std::cerr << "Error: " << e.what() << std::endl;
    return 1;
  }

  return 0;
}
