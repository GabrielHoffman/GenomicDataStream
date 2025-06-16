/***********************************************************************
 * @file  ParallelGenomicChunks.h
 * @author  Zilong Li
 * @brief   get chunks in parallel for pca stream
 * Copyright (C) 2025. The use of this code is governed by the LICENSE file.
 ***********************************************************************/

#ifndef PARALLEL_GENOMIC_CHUNKS_H_
#define PARALLEL_GENOMIC_CHUNKS_H_

#include <atomic>
#include <memory>
#include <type_traits>
#include <algorithm>

#include "GenomicDataStream.h"
#include "GenomicDataStream_virtual.h"
#include "threadpool.hpp"
#include "utils.h"

// Thread-safe chunk processor for coordinating between reader and worker threads
template <typename T>
class MultiGenomicStreamProcessor {
private:
  // Internal class to manage a pool of file readers
  class ReaderPool {
  public:
    ReaderPool(const gds::Param &param, const std::vector<std::string> &regions, size_t numChunks, size_t poolSize) : param(param) {
      // pre-create readers
      for (size_t i = 0; i < poolSize; ++i) {
        readers.push_back(gds::createFileView(param));
      }

      regionSets = gds::chunk_vector(regions, numChunks);
    }

    std::shared_ptr<gds::GenomicDataStream> getReader(size_t i) {
      std::lock_guard<std::mutex> lock(mutex);
      // alter the region so that this reader get a specific chunk
      const std::vector<std::string>& rg = regionSets[i];

      Rcpp::Rcout << "getReader" << std::endl;
      for(auto x: rg){
        Rcpp::Rcout << x << " ";
      }
      Rcpp::Rcout << endl;


      if (readers.empty()) {
        // create a new reader
        auto reader = gds::createFileView(param);
        reader->setRegions(rg);
        return reader;
      }

      auto reader = readers.back();
      readers.pop_back();
      reader->setRegions(rg);
      return reader;
    }

    void returnReader(std::shared_ptr<gds::GenomicDataStream> reader) {
      std::lock_guard<std::mutex> lock(mutex);
      readers.push_back(reader);
    }

  private:
  gds::Param param;
  std::vector<std::vector<std::string> > regionSets;
  std::vector<std::shared_ptr<gds::GenomicDataStream>> readers;
  std::mutex mutex;
  };

public:
  MultiGenomicStreamProcessor(const gds::Param &param, const std::vector<std::string> &regions, size_t numChunks, size_t numThreads)
      : chunkSize(param.chunkSize), 
        numChunks(numChunks),
        pool(numThreads),
        readerPool(param, regions, numChunks, std::max(numThreads, (size_t)std::thread::hardware_concurrency())) // as many as possible readers
  { }

  template<class F>
  void processChunk(F processFunc) {

    std::vector<std::future<void>> futures;

    // submit task to thread pool
    for (size_t i = 0; i < numChunks; ++i) {
      futures.emplace_back(pool.enqueue([this, i, processFunc] {
        // read and process the chunk
        auto reader = readerPool.getReader(i);
        gds::DataChunk<T> chunk;
        try {
          if (reader->getNextChunk(chunk)) {
            processFunc(chunk, i);
          }
          // put the reader back the pool so that we can reuse it later
          readerPool.returnReader(reader);
        } catch (...) {
          // ensure reader is returned to pool even if exception occurs
          readerPool.returnReader(reader);
          Rcpp::Rcout << "reader fails" << endl;
          throw;
        }
      }));
    }

    // wait for all to complete
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



#endif
