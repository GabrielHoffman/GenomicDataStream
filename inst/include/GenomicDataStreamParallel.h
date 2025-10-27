
/*************************************************************
 * @file        GenomicDataStreamParallel.h
 * @author      Gabriel Hoffman
 * @email       gabriel.hoffman@mssm.edu
 * @brief       Read chunks of GenomicStream in parallel
 * Copyright (C) 2025 Gabriel Hoffman
 ************************************************************/

#ifndef MULTI_STREAM_H
#define MULTI_STREAM_H

#if RCPP_PARALLEL_USE_TBB
// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
#else
#include <tbb.h>
#endif

using namespace std;
using namespace gds;

template <typename T>
class GenomicDataStreamParallel {

  public:
  
  GenomicDataStreamParallel(const Param &param, const vector<string> &regions, size_t numChunks, size_t numThreads) :
    numChunks(min(numChunks, regions.size())),
    numThreads(numThreads), 
    limited_arena(numThreads) {


    Rcpp::Rcout << "GenomicDataStreamParallel constructor" << std::endl;

    // reserve space for each reader
    readers.reserve(numThreads);

    // Initialize each reader, run in parallel
    limited_arena.execute([&] {
    tbb::parallel_for(
    tbb::blocked_range<int>(0, numThreads, 1), 
    [&](const tbb::blocked_range<int>& r){ 

        Rcpp::Rcout << "createFileView: " << r.begin() << std::endl;
        readers[r.begin()] = gds::createFileView(param);
        Rcpp::Rcout << "createFileView...done" << std::endl;
      });
    });

    Rcpp::Rcout << "After createFileViews run in parallel" << std::endl;

    // split regions by chunk
    regionSets = gds::chunk_vector(regions, numChunks);
  }

  void processChunks(std::function<void(const gds::DataChunk<T> &, size_t)> processFunc) {

    bool useFilter = true;

    vector<DataChunk<T> > chunkSet;
    chunkSet.reserve(min(numThreads, numChunks));


    Rcpp::Rcout << "numChunks: " << numChunks << std::endl;
    Rcpp::Rcout << "numThreads: " << numThreads << std::endl;
    Rcpp::Rcout << "readers.size(): " << readers.size() << std::endl;


    mutex mtx;  

    // like a thread pool, but just stores index of thread
    // initialize threadSet as 0:numThreads
    vector<int> threadSet( min(numThreads, numChunks));
    iota(threadSet.begin(), threadSet.end(), 0);

    // Parallel part using Thread Building Blocks
    limited_arena.execute([&] {
    tbb::parallel_for(
    tbb::blocked_range<int>(0, numChunks, 1), 
    [&](const tbb::blocked_range<int>& r){ 

        // Get threadIdx from threadSet
        // lock this, so only one thread uses readers[i]
        int threadIdx;
        {
          lock_guard<mutex> lock(mtx);
          threadIdx = threadSet.back();
          threadSet.pop_back();
        }

        int chunkIdx = r.begin(); // which chunk
        Rcpp::Rcout << "chunkIdx: " << chunkIdx << std::endl;
        
        // get reader and set region
        Rcpp::Rcout << "threadIdx: " << threadIdx << std::endl;
        for( auto x: regionSets[chunkIdx]){
          Rcpp::Rcout << x << std::endl;
        }
        auto reader = readers[threadIdx];
        Rcpp::Rcout << "got reader: " << std::endl;
        reader->setRegions( regionSets[chunkIdx] );
        Rcpp::Rcout << "setRegions " << std::endl;

        // Get data and run analysis function
        Rcpp::Rcout << "getchunkSet " << std::endl;
        // DataChunk<T> chunk = chunkSet[threadIdx];

        Rcpp::Rcout << "read variants " << std::endl;
        // while reader get variants
        while (reader->getNextChunk(chunkSet[threadIdx], useFilter)) {
          processFunc(chunkSet[threadIdx], chunkIdx);
        }

        Rcpp::Rcout << "reset mutex " << std::endl;
        // push thread back into the pool
        {
          lock_guard<mutex> lock(mtx);
          threadSet.push_back(threadIdx);
        }

      });
    });
  }

  private:
  int numChunks, numThreads;
  tbb::task_arena limited_arena;
  vector<shared_ptr<GenomicDataStream> > readers;
  vector<vector<string> > regionSets;
};






#endif