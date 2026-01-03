# Changelog

## GenomicDataStream 0.0.66

- add lzma.h

## GenomicDataStream 0.0.65

- update definitions of MatrixInfo and DataChunk to enable filtering

## GenomicDataStream 0.0.64

- add
  [`readH5AD()`](http://gabrielhoffman.github.io/GenomicDataStream/reference/readH5AD.md)
  with `SingleCellExperiment` with file-backed counts

## GenomicDataStream 0.0.63

- Improve performance with no variant filtering

## GenomicDataStream 0.0.62

- Fix bugs in `PCA`, `DataChunk` and `GenomicDataStreamParallel`

## GenomicDataStream 0.0.61

- update Makefile to resolve install issues

## GenomicDataStream 0.0.60

- June 30, 2025
- Faster PCA for VCF and BCF that solves filtering issue

## GenomicDataStream 0.0.59

- June 23, 2025
- `bgenstream` now supports reading a large number of genomic intervals
  using `VariantSet` and then querying by rsid instead of postion,
  instead of crashing at ~10k
- In `VariantSet`, add `id` field
- In `pgenstream` move where `VariantSet` is initialized

## GenomicDataStream 0.0.58

- June 18, 2025
- Article addressing `RcppParallel` on Mac OS X
- Fix Makefile

## GenomicDataStream 0.0.57

- June 16, 2025
- Fix issue with multiple variants at the same site
- PCA is now parallelized with
  [oneTBB](https://uxlfoundation.github.io/oneTBB/)
- add `GenomicDataStreamParallel.h`

## GenomicDataStream 0.0.54

- June 10, 2025
- Move `Depends` to `Imports`

## GenomicDataStream 0.0.53

- June 3, 2025
- [`PCAstream()`](http://gabrielhoffman.github.io/GenomicDataStream/reference/PCAstream.md)
  works on `DelayedArray` and `SingleCellExperiment`
- Add `scaleAndCenter` argument
- add `compute_center_and_scale()`

## GenomicDataStream 0.0.52

- May 30, 2025
- Update
  [`PCAstream()`](http://gabrielhoffman.github.io/GenomicDataStream/reference/PCAstream.md)
- update C++ streams to add `useFilter`
- in C++, add `applyVariantFilter()` based on MAF for SNPs
- Update `DelayedStream`
- PCAOne works for `DelayedStream`

## GenomicDataStream 0.0.51

- May 21, 2025
- Fix bug reading plink files

## GenomicDataStream 0.0.50

- May 19, 2025
- Fix issue with buffer overflow when reading plink files

## GenomicDataStream 0.0.40

- May 15, 2025
- Progress on
  [`PCAstream()`](http://gabrielhoffman.github.io/GenomicDataStream/reference/PCAstream.md)
- parsing plink BIM and FAM files supports tab and space delimiters

## GenomicDataStream 0.0.30

- May 12, 2025
- Add new
  [`PCAstream()`](http://gabrielhoffman.github.io/GenomicDataStream/reference/PCAstream.md)
  using `stream_pcaone()` algorithm with C++ code contributed by
  [@Zilong-Li](https://github.com/Zilong-Li)

## GenomicDataStream 0.0.20

- April 21, 2025
- Resolve `R CMD check` Notes
- refactor for compatibility across ecosystem

## GenomicDataStream 0.0.19

- April 4, 2025
- use `GabrielHoffman/vcfppR` in `DESCRIPTION`

## GenomicDataStream 0.0.18

- April 4, 2025
- bug fix in `main`

## GenomicDataStream 0.0.17

- March 12, 2025
- add `winSVDstream()`

## GenomicDataStream 0.0.16

- Feb 12, 2025
- add `MAF` / `minVariance` filter to
  [`GenomicDataStream()`](http://gabrielhoffman.github.io/GenomicDataStream/reference/GenomicDataStream.md)
  and backend C++ code
- add colnames to `coef` and `se`
- add
  [`getSampleNames()`](http://gabrielhoffman.github.io/GenomicDataStream/reference/getSampleNames.md)

## GenomicDataStream 0.0.15

- Feb 5, 2025
- `GenomicDataStream` can parse genetics file one time, then use
  [`setRegion()`](http://gabrielhoffman.github.io/GenomicDataStream/reference/setRegion.md)
  multiple times to query different regions
- changes to `vcfstream.h`, `pgenstream.h`, `bgenstream.h` to allow
  single initialization and multiple queries with
  [`setRegion()`](http://gabrielhoffman.github.io/GenomicDataStream/reference/setRegion.md)
  multiple times
- fix issue with `vcfstream::setRegions()` failing when empty region was
  used after a valid region

## GenomicDataStream 0.0.14

- Jan 30, 2025
- update `Makevars.in`
- get BH v1.84 from `eddelbuettel/bh@852b468`

## GenomicDataStream 0.0.13

- Jan 29, 2025
- remove dependency on `Rhtslib`

## GenomicDataStream 0.0.10

- Jan 23, 2025
- update to `pgenlibr` v0.4.0
- move regression to separate package
- remove `src/zstd-1.1.0`

## GenomicDataStream 0.0.9

- Dec 18, 2024
- use `fitLinReg` 0.0.5
- check BH (== 1.84.0.0) since 1.87.0-1 is not compatible with boost in
  bgen library

## GenomicDataStream 0.0.8

- Dec 16, 2024
- bug fixes and convergence criteria

## GenomicDataStream 0.0.6

- Dec 3, 2024
- add dependency `fastLinReg`

## GenomicDataStream 0.0.5

- Nov 19, 2024
- Improve interface, documentation and integration with `fastlmm`

## GenomicDataStream 0.0.4

- Oct 31, 2024
- towards support for `pgenlibr`

## GenomicDataStream 0.0.3

- Oct 29, 2024
- `DelayedStream` fixed
- add genetic `Datainfo`

## GenomicDataStream 0.0.2

- Oct 24, 2024
- Add support for BGEN 1.x with phased or unphased biallelic variants
- test data in `inst/extdata/`
- extensive testing for multiple file formats
- Fix BGEN issue with NaN for unphased data

## GenomicDataStream 0.0.1

- Sept 29, 2024
- Add support for multiple file types
