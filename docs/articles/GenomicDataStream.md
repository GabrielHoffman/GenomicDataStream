# Reading data as a GenomicDataStream

`GenomicDataStream` is designed to read chunks of *features* rather than
chunks of *samples*. Features are stored as columns in the matrix
returned by R/C++, *independent of the underlying data storage format*.

## Usage

`GenomicDataStreamRegression` implements regression models (linear and
GLMs) that stream chucks of features using the `GenomicDataStream`
interface. In general, variants from genetic data are used as covariates
in `lmFitFeatures()`, and genes from single cell data are used as
responses in
[`lmFitResponses()`](https://rdrr.io/pkg/BatchRegression/man/lmFitResponses.html).

#### Example code with R

##### Read genotype data into R

``` r

library(GenomicDataStream)
library(GenomicDataStreamRegression) 

# VCF file
file <- system.file("extdata", "test.vcf.gz", package = "GenomicDataStream")

# initialize 
gds = GenomicDataStream(file, "DS", chunkSize=5, initialize=TRUE)

n = 60
y = rnorm(n)
names(y) = paste0("I", seq(n))
design = model.matrix(~1, data.frame(y))

# loop until break
while( 1 ){

  # get data chunk
  # data$X matrix with features as columns
  # data$info information about each feature as rows
  dat = getNextChunk(gds)

  # check if end of stream 
  if( atEndOfStream(gds) ) break

  # do analysis on this chunk of data
  # from GenomicDataStreamRegression
  fit = lmFitFeatures(y, design, dat$X)
}
```

##### Use R to run analysis at C++ level

``` r

library(GenomicDataStream)
library(GenomicDataStreamRegression) 

# VCF file
file <- system.file("extdata", "test.vcf.gz", package = "GenomicDataStream")

# create object, but don't read yet 
# Read DS field storing dosage
gds = GenomicDataStream(file, "DS", chunkSize=5)

n = 60
y = rnorm(n)
names(y) = paste0("I", seq(n))
design = model.matrix(~1, data.frame(y))

# regression of y ~ design + X[,j]
#   where X[,j] is the jth variant in the GenomicDataStream
# data in GenomicDataStream is only accessed at C++ level 
# from GenomicDataStreamRegression
fit = lmFitFeatures(y, design, gds)
```

#### Example code with C++17

``` cpp
#include <RcppArmadillo.h>
#include <GenomicDataStream.h>
 
// use namespace for GenomicDataStream
using namespace gds;
 
// parameters 
string file = "test.vcf.gz";
string field = "DS";    // read dosage field
string region = "";     // no region filter
string samples = "-";   // no samples filter
double MAF = 0.05;   // minor allele freq filter
double minVariance = 0; // retain features with var > minVariance 
int chunkSize = 4;      // each chunk will read 4 variants
 
// initialize parameters
Param param( file, region, samples, MAF, minVariance, chunkSize);
param.setField(field);
 
// Initialise GenomicDataStream to read 
// VCF/BCF/BGEN/PGEN with same interface
shared_ptr<GenomicDataStream> gdsStream = createFileView( param );
 
// declare DataChunk storing an Armadillo matrix for each chunk
DataChunk<arma::mat> chunk;
 
// Store meta-data about each variant
VariantInfo *info;
 
// loop through chunks
while( gdsStream->getNextChunk( chunk ) ){

  // get data from chunk
  arma::mat X = chunk.getData();

  // get variant information
  info = chunk.getInfo<VariantInfo>();

  // Do analysis with variants in this chunk
  analysis_function(X, info);
}
```

Accessing data

## Data types

At the C++ level data from the `GenomicDataStream` can be accessed as
multiple types a according to the templated definition of `DataChunk`.
Above, each `chunk` is a `DataChunk<arma::mat>` so the data is returned
as an Armadillo matrix of doubles (i.e. `arma::mat`).

Suppported data types for dense matricies are:

| Library   | Namespace | Type             |
|-----------|-----------|------------------|
| Rcpp      | `Rcpp`    | `NumericMatrix`  |
| Armadillo | `arma`    | `mat`            |
| Eigen     | `Eigen`   | `MatrixXd`       |
| STL       | `std`     | `vector<double>` |

Importantly, each of these types wraps an array of `double` storing data
in column-major order. The types just provide different interfaces to
the underlying array for use in downstream analysis. In each case, the
raw data is read as an `double` array, a constructor is called for the
requested data type using the array as the underlying data and an object
is returned to the user. In fact, the constructors to the dense matrix
types for Eigen and Armadillo return objects that point to the original
`double` array, without allocating new memory.

`GenomicDataStream` also supports the following sparse matrix types:

| Library   | Namespace | Type           |
|-----------|-----------|----------------|
| Armadillo | `arma`    | `sp_mat`       |
| Eigen     | `Eigen`   | `SparseMatrix` |

In each case, a constructor is called for the requested data type the
converter `double` array to a sparse matrix.

## Session info

    ## R version 4.5.1 (2025-06-13)
    ## Platform: aarch64-apple-darwin23.6.0
    ## Running under: macOS Sonoma 14.7.1
    ## 
    ## Matrix products: default
    ## BLAS/LAPACK: /opt/homebrew/Cellar/openblas/0.3.30/lib/libopenblasp-r0.3.30.dylib;  LAPACK version 3.12.0
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## time zone: America/New_York
    ## tzcode source: internal
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] digest_0.6.39     desc_1.4.3        R6_2.6.1          fastmap_1.2.0     xfun_0.55        
    ##  [6] cachem_1.1.0      knitr_1.51        htmltools_0.5.9   rmarkdown_2.30    lifecycle_1.0.4  
    ## [11] cli_3.6.5         sass_0.4.10       pkgdown_2.2.0     textshaping_1.0.4 jquerylib_0.1.4  
    ## [16] systemfonts_1.3.1 compiler_4.5.1    tools_4.5.1       ragg_1.5.0        bslib_0.9.0      
    ## [21] evaluate_1.0.5    yaml_2.3.12       otel_0.2.0        jsonlite_2.0.0    rlang_1.1.6      
    ## [26] fs_1.6.6          htmlwidgets_1.6.4

\<\>
