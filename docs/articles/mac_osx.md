# Issues on Mac OS X

`GenomicDataStream` runs analyses across multiple cores in parallel
using Threading Building Blocks in
[oneTBB](https://uxlfoundation.github.io/oneTBB/). This is available in
R using [RcppParallel](https://cran.r-project.org/package=RcppParallel).
This package usually installs without an issue, but if it fails you can
use the package after modifying the paths:

##### Install TBB using [brew](https://brew.sh)

``` sh
brew install tbb
```

##### Install RcppParallel with TBB path

``` r

DIR = system("echo $(brew --prefix tbb)", intern=TRUE)
TBB_INC = paste(DIR, "include/", sep='/')
TBB_LIB = paste(DIR, "lib/", sep='/')

.Internal(Sys.setenv("TBB", DIR))
.Internal(Sys.setenv("TBB_INC", TBB_INC))
.Internal(Sys.setenv("TBB_LIB", TBB_LIB))

install.packages("RcppParallel")
```

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
