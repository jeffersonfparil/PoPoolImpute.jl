# PoPoolImpute.jl
Imputation of allele frequency information from pool sequencing genotype data.

|**Laboratory**|**Build Status**|**License**|
|:---:|:---:|:---:|
| <a href="https://adaptive-evolution.biosciences.unimelb.edu.au/"><img src="https://adaptive-evolution.biosciences.unimelb.edu.au/Adaptive%20Evolution%20Logo%20mod.png" width="150"></a> | <a href="https://github.com/jeffersonfparil/PoPoolImpute.jl/actions"><img src="https://github.com/jeffersonfparil/PoPoolImpute.jl/actions/workflows/julia.yml/badge.svg"></a> | [![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0) |

# Usage
`PopPoolImpute.impute(pileup_with_missing::String; window_size::Int=100, model::String=["Mean", "OLS", "RR", "LASSO", "GLMNET"][2], distance::Bool=true, syncx_imputed::String="", threads::Int=2, lines_per_chunk::Int=10_000)::String`

# Inputs
1. pileup_with_missing \[String\]: filename of input pileup file
2. window_size \[Int; default=100\]: number of loci per window
3. model \[String; default="OLS"\]: imputation model to use. Choose from "Mean" (average allele count), "OLS" (ordinary least squares), "RR" (ridge regression), "LASSO" (least absolute shrinkage and selection operator regression), "GLMNET" (elastic-net regression at α=0.5)
4. distance \[Bool; default=true\]: use the first 3 principal components of the pairwise loci distance matrix as an additional covariate
5. syncx_imputed \[String; default=\${pileup_with_missing%.pileup*}-IMPUTED.syncx\]: filename of the imputation output file
6. threads \[Int; default=1\]: number of computing threads or cores to use
7. lines_per_chunk \[Int; default=10,000\]: number of loci per input file of a parallel imputation process

# Output
**Syncx** format (after popoolation2's sync or synchronised pileup file format):
- Column 1:   chromosome or scaffold name
- Column 2:   locus position repeated 7 times corresponding to alleles "A", "T", "C", "G", "INS", "DEL", "N", where "INS" is insertion, "DEL" is deletion, and "N" is unclassified
- Column 3-n: allele counts one column for each pool or population

# Examples

```
# Single-threaded execution
using PoPoolImpute
PoPoolImpute.impute("test.pileup")

# Multi-threaded execution
using Distributed
int_thread_count = length(Sys.cpu_info())-1
Distributed.addprocs(int_thread_count)
@everywhere using PoPoolImpute
PoPoolImpute.impute("test.pileup", window_size=20, threads=2, lines_per_chunk=30)
```

# Details
Performs a simple least squares linear regression to predict missing allele counts per window for each pool with at least one locus with missing data.
- For each pool with missing data we estimate β̂ as:
```
          yₚ = Xₚβ
        → β̂ = inverse(XₚᵀXₚ) (Xₚᵀyₚ).
```

- For each pool with missing data we, imputation is achieved by predicting the missing allele counts:
```
          ŷₘ = XₘB̂.
```
- Where:
    + yₚ is the vector of allele counts of one of the pools with missing data at the loci without missing data (length: mₚ non-missing loci × 7 alleles);
    + Xₚ is the matrix of allele counts of pools without missing data at the loci without missing data in the other pools (dimensions: mₚ non-missing loci × 7 alleles, nₚ pools without missing loci);
    + β̂ is the vector of estimates of the effects of each pool without missing data on the allele counts of one of the pools with missing data (length: nₚ pools without missing loci);
    + inverse() is the Moore-Penrose pseudoinverse if the automatic Julia solver fails;
    + ŷₘ is the vector of imputed allele counts of one of the pools with missing data (length: mₘ missing loci × 7 alleles); and
    + Xₘ is the matrix of allele counts of pools without missing data at the loci with missing data in the other pools (dimensions: mₘ non-missing loci × 7 alleles, nₚ pools without missing loci).

- The imputed allele counts are averaged across the windows sliding one locus at a time.
