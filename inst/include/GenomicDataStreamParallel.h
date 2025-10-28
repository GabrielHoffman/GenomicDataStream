
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

    // reserve space for each reader
    readers.reserve(numThreads);

    mutex mtx;  

    // Initialize each reader, run in parallel
    //  run createFileView() in parallel, but
    //  push to readers in serial with mutex
    limited_arena.execute([&] {
    tbb::parallel_for(
    tbb::blocked_range<int>(0, numThreads, 1), 
    [&](const tbb::blocked_range<int>& r){ 
        // readers[r.begin()] = gds::createFileView(param);

        auto x = createFileView( param );
        {
          lock_guard<mutex> lock(mtx);
          readers.push_back( x );
        }
      });
    });

    // split regions by chunk
    regionSets = gds::chunk_vector(regions, numChunks);
  }

  void processChunks(std::function<void(const gds::DataChunk<T> &, size_t)> processFunc) {

    bool useFilter = true;

    // vector<DataChunk<T> > chunkSet;
    // chunkSet.reserve(min(numThreads, numChunks));

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
        
        // get reader and set region
        auto reader = readers[threadIdx];
        reader->setRegions( regionSets[chunkIdx] );

        // Get data and run analysis function
        DataChunk<T> chk; 
        // while reader get variants
        while (reader->getNextChunk(chk, useFilter)) {
          processFunc(chk, chunkIdx);
        }

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