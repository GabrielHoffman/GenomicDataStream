/***********************************************************************
 * @file		GenomicDataStream.h
 * @author		Gabriel Hoffman
 * @email		gabriel.hoffman@mssm.edu
 * @brief		GenomicDataStream defines an interface to read chunks of data into memory as a matrix.  Supports VCF/VCFGZ/BCF, BGEN, and DelayedArray.  Importing this header gives access to the entire library and gds namespace.
 * Copyright (C) 2024 Gabriel Hoffman
 ***********************************************************************/

/*! \mainpage A scalable interface between data and analysis underneath R
 *
 \image html "GenomicDataStream.png"
<div style="text-align: justify">
Reading genomic data files ([VCF](https://www.ebi.ac.uk/training/online/courses/human-genetic-variation-introduction/variant-identification-and-analysis/understanding-vcf-format/), [BCF](https://samtools.github.io/bcftools/howtos/index.html), [BGEN](https://www.chg.ox.ac.uk/~gav/bgen_format/index.html), [H5AD](https://anndata.readthedocs.io/en/latest/index.html), [DelayedArray](https://bioconductor.org/packages/DelayedArray)) into R/Rcpp in chunks for analysis with [Armadillo](https://doi.org/10.21105/joss.00026) / [Eigen](eigen.tuxfamily.org) / [Rcpp](https://www.rcpp.org) libraries.  Mondern datasets are often too big to fit into memory, and many analyses operate a small chunk features at a time.  Yet in practice, many implementations require the whole dataset stored in memory.  Others pair an analysis with a specific data format (i.e. regresson analysis paired with genotype data from a VCF) in way that the two components can't be separated for use in other applications. 
 </div>
 The `GenomicDataStream` C++ interface separates 

 -# data source 
 -# streaming chunks of features into a data matrix
 -# downstream analysis  
 * 
 * ### Example code with C++17
```cpp
#include <RcppArmadillo.h>
#include <GenomicDataStream.h>

// use namespace for GenomicDataStream
using namespace gds;

// parameters 
string file = "test.vcf.gz";
string field = "DS";    // read dosage field
string region = "";     // no region filter
string samples = "-";   // no samples filter
int chunkSize = 4;      // each chunk will read 4 variants

// initialize parameters
Param param(file, region, samples, chunkSize);
param.setField( field );

// Initialise GenomicDataStream to read 
// VCF/BCF and BGEN with same interface
unique_ptr<GenomicDataStream> gdsStream = createFileView( param );

// declare DataChunk storing an Armadillo matrix for each chunk
DataChunk<arma::mat> chunk;

// Store meta-data about each variant
VariantInfo *info;

// loop through chunks
while( gdsStream->getNextChunk( chunk ) ){

    // get data from chunk
    // chunk.getData();

    // get variant information
    info = chunk.getInfo<VariantInfo>();

    // Do analysis with variants in this chunk
}
```

 
 * ## Dependencies
| Package | Ref | Role |
| - | --- | --------- |
[vcfppR](https://cran.r-project.org/package=vcfppR) | [Bioinformatics](https://doi.org/10.1093/bioinformatics/btae049)  | C++ API for htslib  |
[htslib](https://github.com/samtools/htslib) | [GigaScience](https://doi.org/10.1093/gigascience/giab007)  | C API for VCF/BCF files |
[beatchmat](https://bioconductor.org/packages/beachmat/) | [PLoS Comp Biol](https://doi.org/10.1371/journal.pcbi.1006135)  | C++ API for access data owned by R |
[Rcpp](https://cran.r-project.org/package=Rcpp)| [J Stat Software](https://doi.org/10.18637/jss.v040.i08) |  API for R/C++ integration
[RcppEigen](https://cran.r-project.org/package=RcppEigen) | [J Stat Software](https://doi.org/10.18637/jss.v052.i05) | API for Rcpp access to Eigen matrix library
[RcppArmadillo](https://cran.r-project.org/package=RcppArmadillo)| [J Stat Software](https://doi.org/10.18637/jss.v040.i08) | API for Rcpp access to Armadillo matrix library
[Eigen](https://eigen.tuxfamily.org) | |C++ library for linear algebra with advanced features
[Armadillo](https://arma.sourceforge.net) | [J Open Src Soft](https://doi.org/10.21105/joss.00026) | User-friendly C++ library for linear algebra
 * 
 * 
## Notes
<div style="text-align: justify">
`GenomicDataStream` provide flexability in terms of data input types and and matrix libraries.  This can useful in many cases, but the large number of dependencies can require installation of additional libraries and increase compile times.  Some of these dependencies can be avoided by removing support for some capabilities with compiler flags in `Makevars`:
</div>

 `-D DISABLE_DELAYED_STREAM`     
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Omit `DelayedStream` class, remove dependence on `Rcpp` and `beachmat`  
 
 `-D DISABLE_EIGEN`   
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Omit support for Eigen matrix library, and remove dependence on `RcppEigen` and `Eigen`

 `-D DISABLE_RCPP`  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 
Omit support for `Rcpp` matrix library, and remove dependence on `Rcpp`
 */


#ifndef GENOMIC_DATA_STREAM_H_
#define GENOMIC_DATA_STREAM_H_

#include "pgenstream.h"
#include "bgenstream.h"
#include "vcfstream.h"
#include "DelayedStream.h"
#include "utils.h"
#include "MatrixInfo.h"
#include "VariantInfo.h"

namespace gds {

/** Create a reader for VCF/VCFGZ/BCF and BGEN that instantiates either vcfstream or bgenstream based on the file extension in param.
 * @param param class of Param type
 */ 
static unique_ptr<GenomicDataStream> createFileView( const Param & param ){

    unique_ptr<GenomicDataStream> gdsStream;

    // Define reader for VCF/VCFGZ/BCF or BGEN
    // depending on file extension
    switch( param.fileType ){
        case VCF:
        case VCFGZ:
        case BCF:
            gdsStream = make_unique<vcfstream>( param );
            break;
        case BGEN:
            gdsStream = make_unique<bgenstream>( param );
            break;
        case PGEN:
        case PBED:
            gdsStream = make_unique<pgenstream>( param );
            break;
        case OTHER:
            throw runtime_error("Invalid file extension: " + param.file);
            break;
    }  

    return gdsStream; 
}

/** Create a reader for VCF/VCFGZ/BCF and BGEN that instantiates either vcfstream or bgenstream based on the file extension in param.
 * @param param class of Param type
 */ 
static shared_ptr<GenomicDataStream> createFileView_shared( const Param & param ){

    shared_ptr<GenomicDataStream> gdsStream;

    // Define reader for VCF/VCFGZ/BCF or BGEN
    // depending on file extension
    switch( param.fileType ){
        case VCF:
        case VCFGZ:
        case BCF:
            gdsStream = make_shared<vcfstream>( param );
            break;
        case BGEN:
            gdsStream = make_shared<bgenstream>( param );
            break;
        case PGEN:
        case PBED:
            gdsStream = make_unique<pgenstream>( param );
            break;
        case OTHER:
            throw runtime_error("Invalid file extension: " + param.file);
            break;
    }  

    return gdsStream; 
}

/* Defines type for interface with Rcpp */
typedef struct BoundDataStream {
    BoundDataStream(const Param &param){
        ptr = createFileView_shared( param );
    }

    shared_ptr<gds::GenomicDataStream> ptr;
    bool atEndOfStream = false;
    long featuresRead = 0;
} BoundDataStream;


}

#endif
