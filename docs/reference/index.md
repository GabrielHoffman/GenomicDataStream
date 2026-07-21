# Package index

## GenomicDataStream

- [`GenomicDataStream()`](http://gabrielhoffman.github.io/GenomicDataStream/reference/GenomicDataStream.md)
  : GenomicDataStream
- [`GenomicDataStream-class`](http://gabrielhoffman.github.io/GenomicDataStream/reference/GenomicDataStream-class.md)
  : Interface to genomic data files
- [`getNextChunk()`](http://gabrielhoffman.github.io/GenomicDataStream/reference/getNextChunk.md)
  : Get data chunk from GenomicDataStream
- [`initializeStream()`](http://gabrielhoffman.github.io/GenomicDataStream/reference/initializeStream.md)
  : Initialize GenomicDataStream
- [`setChunkSize()`](http://gabrielhoffman.github.io/GenomicDataStream/reference/setChunkSize.md)
  : Set Chunk Size
- [`setRegion()`](http://gabrielhoffman.github.io/GenomicDataStream/reference/setRegion.md)
  : Set regions of GenomicDataStream

## Read H5AD file

- [`readH5AD()`](http://gabrielhoffman.github.io/GenomicDataStream/reference/readH5AD.md)
  : Read H5AD as SingleCellExperiment

## Out-of-core SVD

- [`dashSVD()`](http://gabrielhoffman.github.io/GenomicDataStream/reference/dashSVD.md)
  : Dynamic Shifts Based Randomized SVD
- [`PCAone()`](http://gabrielhoffman.github.io/GenomicDataStream/reference/PCAone.md)
  : Window-based Randomized SVD with PCAone
- [`PCA-class`](http://gabrielhoffman.github.io/GenomicDataStream/reference/PCA-class.md)
  : PCA result

## Accessor functions

- [`atEndOfStream()`](http://gabrielhoffman.github.io/GenomicDataStream/reference/atEndOfStream.md)
  : Detected if end of stream is reached
- [`countChunks()`](http://gabrielhoffman.github.io/GenomicDataStream/reference/countChunks.md)
  : Get chunk sizes
- [`featuresRead()`](http://gabrielhoffman.github.io/GenomicDataStream/reference/featuresRead.md)
  : Get number of features read from GenomicDataStream
- [`getSampleNames()`](http://gabrielhoffman.github.io/GenomicDataStream/reference/getSampleNames.md)
  : Get sample names
- [`getVariantLocations()`](http://gabrielhoffman.github.io/GenomicDataStream/reference/getVariantLocations.md)
  : Get location of each feature
- [`isInitialized()`](http://gabrielhoffman.github.io/GenomicDataStream/reference/isInitialized.md)
  : Get status of GenomicDataStream
- [`plot(`*`<PCA>`*`,`*`<ANY>`*`)`](http://gabrielhoffman.github.io/GenomicDataStream/reference/plot-methods.md)
  : Plot PCAone
- [`print(`*`<GenomicDataStream>`*`)`](http://gabrielhoffman.github.io/GenomicDataStream/reference/print-methods.md)
  [`print(`*`<PCA>`*`)`](http://gabrielhoffman.github.io/GenomicDataStream/reference/print-methods.md)
  : Print object
- [`rownames(`*`<GenomicDataStream>`*`)`](http://gabrielhoffman.github.io/GenomicDataStream/reference/rownames.md)
  : Get rownames
- [`show(`*`<GenomicDataStream>`*`)`](http://gabrielhoffman.github.io/GenomicDataStream/reference/show-methods.md)
  [`show(`*`<PCA>`*`)`](http://gabrielhoffman.github.io/GenomicDataStream/reference/show-methods.md)
  : Show object
- [`getChromRanges()`](http://gabrielhoffman.github.io/GenomicDataStream/reference/getChromRanges.md)
  : Get ranges for each chromosome
- [`chopChroms()`](http://gabrielhoffman.github.io/GenomicDataStream/reference/chopChroms.md)
  : Get genomic intervals
- [`summarizeChunks()`](http://gabrielhoffman.github.io/GenomicDataStream/reference/summarizeChunks.md)
  : Get information about chunks

## Other functions

- [`standardize_in_place()`](http://gabrielhoffman.github.io/GenomicDataStream/reference/standardize_in_place.md)
  : Standardize matrix columns in place
- [`normPC()`](http://gabrielhoffman.github.io/GenomicDataStream/reference/normPC.md)
  : Normalize principal components
- [`perfMetric()`](http://gabrielhoffman.github.io/GenomicDataStream/reference/perfMetric.md)
  : Evaluate performance of PC estimates
