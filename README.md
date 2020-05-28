# HAPPI GWAS

<!-- badges: start -->
<!-- badges: end -->

HAPPI GWAS is a genome-wide association study (GWAS) tool written in R.

## Installation

You can install the HAPPI.GWAS from [Github](https://github.com/Angelovici-Lab/HAPPI.GWAS) with:

``` r
# Run this inside R environment
install.packages("devtools", dependencies = TRUE)
devtools::install_github("Angelovici-Lab/HAPPI.GWAS")
```

``` r
# Run this in your terminal
mkdir HAPPI_GWAS
cd HAPPI_GWAS
git clone https://github.com/Angelovici-Lab/HAPPI.GWAS.git
cd HAPPI.GWAS
Rscript setup.R
```

## Usage

``` r
Rscript HAPPI_GWAS.R [-h] [-generateBLUP] [-generateBLUE] [-GAPIT] [-extractHaplotype] [-searchGenes] input

positional arguments:
  input              Input YAML File

optional arguments:
  -h, --help         show this help message and exit
  -generateBLUP      Generate BLUP data from raw data
  -generateBLUE      Generate BLUE data from raw data
  -GAPIT             Run GAPIT
  -extractHaplotype  Extract haplotype (Require: -GAPIT)
  -searchGenes       Search genes (Require: -GAPIT)
```

## Example

This is a basic example which shows you how to use HAPPI.GWAS:

``` r
cd /path/to/HAPPI_GWAS/HAPPI.GWAS
Rscript HAPPI_GWAS.R -GAPIT -extractHaplotype -searchGenes Demo_GLM.yaml
Rscript HAPPI_GWAS.R -GAPIT -extractHaplotype -searchGenes Demo_FarmCPU.yaml
```

## Package Update

To upgrade HAPPI.GWAS to the latest version, please remove the package and re-install the latest HAPPI.GWAS package:

``` r
# Run this inside R environment
remove.packages("HAPPI.GWAS")
devtools::install_github("Angelovici-Lab/HAPPI.GWAS")
```

``` r
# Run this in your terminal
cd /path/to/HAPPI_GWAS/HAPPI.GWAS
git pull
```

## Wiki

[HAPPI GWAS Wiki Page](https://github.com/Angelovici-Lab/HAPPI.GWAS/wiki)

## User Manual

[HAPPI GWAS User Manual](https://github.com/Angelovici-Lab/HAPPI.GWAS/wiki/files/user_manual/HAPPI_GWAS_manual.pdf)

## Updates

05-28-2020: Added HAPPI GWAS user manual   
05-16-2020: Added HAPPI GWAS wiki page 
