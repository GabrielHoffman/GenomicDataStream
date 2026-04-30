# Extending GenomicDataStream

## Additional file types

Consider extending the package to support to support a new file of type
`X`. First, create a new class `Xstream` that inherits from
`GenomicDataStream` defined
[here](http://gabrielhoffman.github.io/GenomicDataStream/doxygen/html/classgds_1_1_genomic_data_stream.md).
`Xstream` must define a constructor taking `Param` defined
[here](http://gabrielhoffman.github.io/GenomicDataStream/doxygen/html/structgds_1_1_param.md),
and implement interfaces for
[`getNextChunk()`](http://gabrielhoffman.github.io/GenomicDataStream/reference/getNextChunk.md)
for each data type.

A
[`DataChunk`](http://gabrielhoffman.github.io/GenomicDataStream/doxygen/html/classgds_1_1_data_chunk.md)
stores information about each feature and sample in a class inheriting
from
[`DataInfo`](http://gabrielhoffman.github.io/GenomicDataStream/doxygen/html/classgds_1_1_data_info.md).
Streams for genotype data store feature information in the
[`VariantInfo`](http://gabrielhoffman.github.io/GenomicDataStream/doxygen/html/classgds_1_1_variant_info.md)
class, while `DelayedStream` uses
[`MatrixInfo`](http://gabrielhoffman.github.io/GenomicDataStream/doxygen/html/classgds_1_1_matrix_info.md).
These classes or a custom class can be used to store information for a
new file type.

While nothing else about the underlying data form at is assumed,
`GenomicDataStream` designed to chunks of *features* rather than chunks
of *samples*. Features are stored as columns in the matrix returned by
R/C++, *independent of the underlying data storage format*.

## Additional matrix types

First, take a look a currently supported [data
types](http://gabrielhoffman.github.io/GenomicDataStream/articles/GenomicDataStream.html#data-types),
and note that each of these types wraps an array of `double` storing
data in column-major order. In fact, the constructors to the dense
matrix types for Eigen and Armadillo return objects that point to the
original `double` array, without allocating new memory. So any new data
type should have a constructor that takes a `double` array.

The simplest new type to implement would be a `float` version of an
existing type. This could reduce memory usage in downstream analyses and
increase speed at the cost of numerical precision.

## Importing into other projects

`GenomicDataStream` is written so that core functions are in C++17 with
no dependence or R or Rcpp. On top of that, there is a thin wrapper that
uses Rcpp to interface between R and the lower-level library.

The C++ code is divided into two sections:

- `inst/include`: header-only C++17 code with no R dependencies

- `src/export.cpp`: `Rcpp` layer between header-only library to R
  interface

`GenomicDataStream` provides flexability in terms of data input types
and and matrix libraries. This can useful in many cases, but the large
number of dependencies can require installation of additional libraries
and increase compile times. Some of these dependencies can be avoided by
removing support for some capabilities with compiler flags in
`Makevars`:

`-D DISABLE_DELAYED_STREAM`\
           Omit `DelayedStream` class, remove dependence on `Rcpp` and
`beachmat`

`-D DISABLE_EIGEN`\
           Omit support for Eigen matrix library, and remove dependence
on `RcppEigen` and `Eigen`

`-D DISABLE_RCPP`\
           Omit support for `Rcpp` matrix library, and remove dependence
on `Rcpp`

`-D DISABLE_PLINK`\
           Omit support for `PLINK` files (PGEN, BED), and remove
dependence on `pgenlibr`
