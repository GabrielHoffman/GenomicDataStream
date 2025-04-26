/***********************************************************************
 * @file		ParallelGenomicChunks.h
 * @author		Zilong Li
 * @brief		get chunks in parallel for pca stream
 * Copyright (C) 2025. The use of this code is governed by the LICENSE file.
 ***********************************************************************/

#include <atomic>
#include <fstream>
#include <iostream>
#include <memory>
#include <type_traits>
#include <vector>

#include <RcppEigen.h>

#include "GenomicDataStream.h"
#include "GenomicDataStream_virtual.h"
#include "threadpool.hpp"

// Thread-safe chunk cache for coordinating between reader and worker threads
template <typename T>
class MultiGenomicStreamProcessor {
private:
  // Internal class to manage a pool of file readers
  class ReaderPool {
  public:
    ReaderPool(const gds::Param &param, size_t poolSize) : param(param) {
      // Pre-create readers
      for (size_t i = 0; i < poolSize; ++i) {
        readers.push_back(gds::createFileView_shared(param));
        // Don't open them yet to avoid having too many open file handles
      }
    }

    std::shared_ptr<gds::GenomicDataStream> getReader() {
      std::lock_guard<std::mutex> lock(mutex);

      if (readers.empty()) {
        // If pool is exhausted, create a new reader
        return gds::createFileView_shared(param);
      }

      auto reader = readers.back();
      readers.pop_back();
      return reader;
    }

    void returnReader(std::shared_ptr<gds::GenomicDataStream> reader) {
      std::lock_guard<std::mutex> lock(mutex);
      readers.push_back(reader);
    }

  private:
    gds::Param param;
    std::vector<std::shared_ptr<gds::GenomicDataStream>> readers;
    std::mutex mutex;
  };

public:
  MultiGenomicStreamProcessor(const gds::Param &param, size_t numThreads)
      : chunkSize(param.chunkSize), numChunks(param.regions.size()),
        pool(numThreads),
        readerPool(param, std::min(numThreads, size_t(16))) // max 16 readers
  {

    std::cout << "Processing in " << numChunks << " chunks\n";
  }

  void
  processFile(std::function<void(const gds::DataChunk<T> &, size_t)> processFunc) {

    std::vector<std::future<void>> futures;

    // Submit task to thread pool
    for (size_t i = 0; i < numChunks; ++i) {
      futures.emplace_back(pool.enqueue([this, i, processFunc] {
        // Get a reader from the pool
        auto reader = readerPool.getReader();
        gds::DataChunk<T> chunk;
        try {
          // read and process the chunk
          if (reader->getNextChunk(chunk)) {
            processFunc(chunk, i);
          }
          // Return the reader to the pool
          readerPool.returnReader(reader);
        } catch (...) {
          // Ensure reader is returned to pool even if exception occurs
          readerPool.returnReader(reader);
          throw;
        }
      }));
    }

    // Wait for all processing to complete
    for (auto &future : futures) {
      future.get();
    }
  }

private:
  const size_t chunkSize;
  const size_t numChunks;
  ThreadPool pool;
  ReaderPool readerPool;
};

int get_chunks_parallel_vcf(const std::string &file,
                            const std::string &region) {

  std::string field = "GT";  // read dosage field
  std::string samples = "-"; // no samples filter
  double minVariance = NAN;  // retain features with var > minVariance
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
    MultiGenomicStreamProcessor<Eigen::MatrixXd> processor(param, 4);

    processor.processFile(
        [&](const gds::DataChunk<Eigen::MatrixXd> &chunk, size_t chunkIndex) {
          // Example processing logic
          auto X = chunk.getData();
          size_t size = X.cols();

          // Calculate sum of X for test
          double chunkSum = X.sum();

          // Thread-safe update
          if(chunkIndex % 4 == 0)
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
