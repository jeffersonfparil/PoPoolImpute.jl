# PoPoolImpute.jl
Imputation of allele frequency information from pool sequencing genotype data.

|**Laboratory**|**Build Status**|**License**|
|:---:|:---:|:---:|
| <a href="https://adaptive-evolution.biosciences.unimelb.edu.au/"><img src="https://adaptive-evolution.biosciences.unimelb.edu.au/Adaptive%20Evolution%20Logo%20mod.png" width="150"></a> | <a href="https://github.com/jeffersonfparil/PoPoolImpute.jl/actions"><img src="https://github.com/jeffersonfparil/PoPoolImpute.jl/actions/workflows/julia.yml/badge.svg"></a> | [![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0) |

# Usage
```
PopPoolImpute.impute(str_filename_input; n_int_window_size=10, n_flt_maximum_fraction_of_pools_with_missing=0.5, n_flt_maximum_fraction_of_loci_with_missing=0.5, str_filename_output="output-imputed.syncx", n_int_thread_count=1)
```

# Inputs
1. *str_filename_input* [String]: filename of the genotype data in [pileup format (.pileup)](http://samtools.sourceforge.net/pileup.shtml)
2. *n_int_window_size* [Integer; default=10]: size of the sliding window across which imputation is performed
3. *n_flt_maximum_fraction_of_pools_with_missing* [Float; default=0.5]: maximum tolerable fraction of the pools with at least one missing locus
4. *n_flt_maximum_fraction_of_loci_with_missing* [Float; default=0.5]: maximum tolerable fraction of the loci with missing data
5. *bool_use_distance_matrix* [Boolean; default=false]: perform ordinary least squares regression with pairwise loci distances covariate
6. *str_model* [String; default="OLS"]: imputation model to use. Choose from "Mean" (average counts in the non-missing pools), "OLS" (ordinary least squares), "RR" (ridge regression), "LASSO" (least absolute shrinkage and selection operator), and "GLMNET" (elastic net with additional parameter flt_glmnet_alpha; default=0.5)

7. *n_int_thread_count* [Integer; default=1]: number of computing threads to use

# Output
*str_filename_output* [String; default="output-imputed.syncx"; comma-separated file]

**Syncx format** (after popoolation2's sync or synchronised pileup file format):
- *Column 1*:   chromosome or scaffold name
- *Column 2*:   locus position repeated 7 times corresponding to alleles "A", "T", "C", "G", "INS", "DEL", "N", where "INS" is insertion, "DEL" is deletion, and "N" is unclassified
- *Columns 3 to n*: are the allele counts one column for each pool or population

# Installation
```
using Pkg
Pkg.add(url="https://github.com/jeffersonfparil/PoPoolImpute.jl.git")
```

# Examples
```
# Single-threaded execution
using PoPoolImpute
str_filename_input = "test.pileup"
PoPoolImpute.impute(str_filename_input)

# Multi-threaded execution
using Distributed
n_int_thread_count = length(Sys.cpu_info())-1
Distributed.addprocs(n_int_thread_count)
@everywhere using PoPoolImpute
PoPoolImpute.impute(str_filename_input, str_filename_output="multithreaded.syncx", n_int_thread_count=n_int_thread_count)
```
# Details

Performs simple mean extrapolation and various linear regression algorithms to predict missing allele counts per window:

**Xₚₘ = XₚₚB**

**→ B̂ = F(Xₚₚ, Xₚₘ, ...)**.

Imputation is achieved by predicting the missing allele counts:

**X̂ₘₘ = XₘₚB̂**.

Where:

- **Xₚₘ** is the matrix of allele counts of pools with missing data at the loci without missing data (dimensions: **mₚ** non-missing loci × 7 alleles, **nₘ** pools with missing loci);
- **Xₚₚ** is the matrix of allele counts of pools without missing data at the loci without missing data in the other pools (dimensions: **mₚ** non-missing loci × 7 alleles, **nₚ** pools without missing loci);
- **B̂** is the matrix of estimates of the effects of each pool without missing data on the allele counts of the pools with missing data (dimensions: **nₚ** pools without missing loci, **nₘ** pools with missing loci);
- **F()** refers to any of the imputation models implemented, i.e. means of th non-missing pools, ordinary least squares, ridge regression, LASSO, and elastic net regularisation;
- **X̂ₘₘ** is the matrix of imputed allele counts of pools with missing data (dimensions: **mₘ** missing loci × 7 alleles, **nₘ** pools with missing loci); and
- **Xₘₚ** is the matrix of allele counts of pools without missing data at the loci with missing data in the other pools (dimensions: **mₘ** non-missing loci × 7 alleles, **nₚ** pools without missing loci).

The imputed allele counts are averaged across the windows sliding one locus at a time.
