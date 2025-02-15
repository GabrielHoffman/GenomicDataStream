
# GenomicDataStream 0.0.16
 - Feb 12, 2025
 - add `MAF` / `minVariance` filter to `GenomicDataStream()` and backend C++ code
 - add colnames to `coef` and `se`
 - add `getSampleNames()`

# GenomicDataStream 0.0.15
 - Feb 5, 2025
 - `GenomicDataStream` can parse genetics file one time, then use `setRegion()` multiple times to query different regions
  - changes to `vcfstream.h`, `pgenstream.h`, `bgenstream.h` to allow single initialization and multiple queries with `setRegion()` multiple times
 - fix issue with `vcfstream::setRegions()` failing when empty region was used after a valid region


# GenomicDataStream 0.0.14
 - Jan 30, 2025
 - update `Makevars.in`
 - get BH v1.84 from `eddelbuettel/bh@852b468`

# GenomicDataStream 0.0.13
 - Jan 29, 2025
 - remove dependency on `Rhtslib`

# GenomicDataStream 0.0.10
 - Jan 23, 2025
 - update to `pgenlibr` v0.4.0
 - move regression to separate package
 - remove `src/zstd-1.1.0`

# GenomicDataStream 0.0.9
 - Dec 18, 2024
 - use `fitLinReg` 0.0.5
 - check BH (== 1.84.0.0) since 1.87.0-1 is not compatible with boost in bgen library

# GenomicDataStream 0.0.8
 - Dec 16, 2024
 - bug fixes and convergence criteria

# √ 0.0.7
 - Dec 12, 2024
 - add `glmFitResponses()` and `glmFitFeatures()`

# GenomicDataStream 0.0.6
 - Dec 3, 2024
 - add dependency `fastLinReg`

# GenomicDataStream 0.0.5
 - Nov 19, 2024
 - Improve interface, documentation and integration with `fastlmm`


# GenomicDataStream 0.0.4
 - Oct 31, 2024
 - towards support for `pgenlibr`


# GenomicDataStream 0.0.3
 - Oct 29, 2024
 - `DelayedStream` fixed
 - add genetic `Datainfo`


# GenomicDataStream 0.0.2
 - Oct 24, 2024
 - Add support for BGEN 1.x with phased or unphased biallelic variants
 - test data in `inst/extdata/`
 - extensive testing for multiple file formats
  - Fix BGEN issue with NaN for unphased data


# GenomicDataStream 0.0.1
 - Sept 29, 2024
 - Add support for multiple file types

