# HAPPI.GWAS

<!-- badges: start -->
<!-- badges: end -->

A genome-wide association study (GWAS) tool written in R.

## Installation

You can install the HAPPI.GWAS from [Github](https://github.com/Angelovici-Lab/HAPPI.GWAS) with:

``` r
# Run this inside R environment
install.packages("devtools", dependencies = TRUE)
devtools::install_github("Angelovici-Lab/HAPPI.GWAS")
```

``` r
# Run this on your terminal
mkdir HAPPI_GWAS
cd HAPPI_GWAS
git clone https://github.com/Angelovici-Lab/HAPPI.GWAS.git
```

## Example

This is a basic example which shows you how to use HAPPI.GWAS:

``` r
cd HAPPI_GWAS/HAPPI.GWAS
<Edit yaml file>
Rscript setup.R
Rscript HAPPI_GWAS.R -GAPIT -extractHaplotype -searchGenes Demo_GLM.yaml
Rscript HAPPI_GWAS.R -GAPIT -extractHaplotype -searchGenes Demo_FarmCPU.yaml
```
