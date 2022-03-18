module PoPoolImpute

include("functions.jl")
using .functions: SPLIT, IMPUTE
using Dates
using ProgressMeter

### Documentation
"""
# ____________________________________________________________________
# PoPoolImpute: imputation of population- and pool-level genotype data
# Usage
`PopPoolImpute.impute(pileup_with_missing::String; window_size::Int=100, model::String=["Mean", "OLS", "RR", "LASSO", "GLMNET"][2], distance::Bool=true, syncx_imputed::String="", threads::Int=2, lines_per_chunk::Int=10_000)::String`

# Inputs
1. pileup_with_missing [String]: filename of input pileup file
2. window_size [Int; default=100]: number of loci per window which should be around the same size as the sequencing read length (Note: at least 1 pool without any missing data across the window is needed for imputation)
3. model [String; default="OLS"]: imputation model to use. Choose from "Mean" (average allele count), "OLS" (ordinary least squares), "RR" (ridge regression), "LASSO" (least absolute shrinkage and selection operator regression), "GLMNET" (elastic-net regression at α=0.5)
4. distance [Bool; default=true]: use the first 3 principal components of the pairwise loci distance matrix as an additional covariate
5. syncx_imputed [String; default=\${pileup_with_missing%.pileup*}-IMPUTED.syncx]: filename of the imputation output file
6. threads [Int; default=1]: number of computing threads or cores to use
7. lines_per_chunk [Int; default=10,000]: number of loci per input file of a parallel imputation process

# Output
syncx_imputed:
Syncx format (after popoolation2's sync or synchronised pileup file format):
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
Performs simple linear regression to predict missing allele counts per window for each pool with at least one locus with missing data. This imputation method requires at least one pool without missing data across the window. It follows that to maximise the number of loci we can impute, we need to impose a maximum window size equal to the length of the sequencing read used to generate the data, e.g. 100 bp to 150 bp for Illumina reads.

- For each pool with missing data we estimate β̂ as:
```
          yₚ = Xₚβ
        → β̂ = inverse(XₚᵀXₚ) (Xₚᵀyₚ).
```

- For each pool with missing data, imputation is achieved by predicting the missing allele counts:
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

# Author
- Jeff Paril (jeffersonparil@gmail.com; https://orcid.org/0000-0002-5693-4123)
...
"""
function impute(pileup_with_missing::String; window_size::Int=100, model::String=["Mean", "OLS", "RR", "LASSO", "GLMNET"][2], distance::Bool=true, syncx_imputed::String="", threads::Int=1, lines_per_chunk::Int=10_000)::String
    ### Opening remark
    println("")
    println("####################################################################")
    println("PoPoolImpute: imputation of population- and pool-level genotype data")
    println("(version 0.1.0; release 2022/03/10)")
    println("####################################################################")
    println("Input parameters:")
    @show pileup_with_missing
    @show window_size
    @show model
    @show distance
    if threads > 1
        @show threads
        @show lines_per_chunk
    end
    println("####################################################################")
    println(string("Start time: ", Dates.format(now(), "Y-u-dd(E)THH:MM")))
    println("####################################################################")
    ### Define output file if not specified
    if syncx_imputed == ""
        syncx_imputed = string(join(split(pileup_with_missing, '.')[1:(end-1)], '.'), "-IMPUTED.syncx")
    end
    ### Define the full path to the input and output files since calling functions within @distributed loop will revert back to the root directory from where julia was executed from
    if dirname(pileup_with_missing) == ""
        pileup_with_missing = string(pwd(), "/", pileup_with_missing)
    end
    if dirname(syncx_imputed) == ""
        syncx_imputed = string(pwd(), "/", syncx_imputed)
    end
    if threads > 1
        ### Multi-threaded execution
        ### Split input file for parallel processing if we have more than 1 core or thread available
        filenames = [pileup_with_missing]
        if threads > 1
            filenames = SPLIT(pileup_with_missing, lines_per_chunk, window_size)
        end
        ### Impute
        @time filenames_out = @sync @showprogress @distributed (append!) for f in filenames
            filename_imputed = IMPUTE(f, window_size=window_size, model=model, distance=distance)
            [filename_imputed]
        end
        ### Sort the chunks so we maintain the one-to-one correspondence between input and output loci arrangement
        sort!(filenames_out)
        ### Trim-off overhanging windows and merge
        file_out = open(syncx_imputed, "w")
        for i in 1:length(filenames_out)
            if i < length(filenames_out)
                ### trim trailing window from the first and intervening chunks
                lines = 0
                file_in = open(filenames_out[i], "r")
                while !eof(file_in)
                    lines += 1
                    readline(file_in);
                end
                close(file_in)
                max_line = lines - window_size
            else
                ### do not trim the last chunk
                max_line = Inf
            end
            file_in = open(filenames_out[i], "r")
            j = 0
            while (!eof(file_in)) & (j < max_line)
                j += 1
                write(file_out, string(readline(file_in), "\n"))
            end
            close(file_in)
            ### clean up
            rm(filenames[i])        ### pileup chunks
            rm(filenames_out[i])    ### syncx chunks
        end
        close(file_out)
    else
        ### Single-threaded execution
        @show pwd()
        @show pileup_with_missing
        syncx_imputed = IMPUTE(pileup_with_missing, window_size=window_size, model=model, distance=distance)
    end
    ### Closing remark
    println("")
    println("####################################################################")
    println(string("End time: ", Dates.format(now(), "Y-u-dd(E)THH:MM")))
    println("Imputation successful. Please find the output file:")
    println(syncx_imputed)
    println("####################################################################")
    return(syncx_imputed)
end

end
