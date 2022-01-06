module PoPoolImpute

### Load linear algebra library for the Moore-Penrose pseudoinverse if the automatic solver fails
using Distributed
using ProgressMeter
using LinearAlgebra
### Load the functions, and move them into scope
include("functions.jl")
using .functions: fun_ascii_allele_states_to_counts_per_locus,
                                fun_ascii_allele_states_to_counts_per_window, 
                                fun_impute_per_window, 
                                fun_simple_progress_bar,
                                fun_writeout_inrun, 
                                fun_single_threaded_imputation
### Documentation
"""
# ____________________________________________________________________
# PoPoolImpute: imputation of population- and pool-level genotype data


# Usage
`PopPoolImpute.impute(str_filename_input; n_int_window_size=10, n_flt_maximum_fraction_of_pools_with_missing=0.5, n_flt_maximum_fraction_of_loci_with_missing=0.5, str_filename_output="output-imputed.syncx", n_int_thread_count=1)`


# Inputs
1. str\\_filename\\_input [String]: filename of the genotype data in [pileup format (.pileup)](http://samtools.sourceforge.net/pileup.shtml)
2. n\\_int\\_window\\_size [Integer; default=10]: size of the sliding window across which imputation is performed
3. n\\_flt\\_maximum\\_fraction\\_of\\_pools\\_with\\_missing [Float; default=0.5]: maximum tolerable fraction of the pools with at least one missing locus
4. n\\_flt\\_maximum\\_fraction\\_of\\_loci\\_with\\_missing [Float; default=0.5]: maximum tolerable fraction of the loci with missing data
5. n\\_int\\_thread\\_count [Integer; default=1]: number of computing threads to use


# Output
str\\_filename\\_output [String; default="output-imputed.syncx"; comma-separated file]

Syncx format (after popoolation2's sync or synchronised pileup file format):
- Column 1:   chromosome or scaffold name
- Column 2:   locus position repeated 7 times corresponding to alleles "A", "T", "C", "G", "INS", "DEL", "N", where "INS" is insertion, "DEL" is deletion, and "N" is unclassified
- Column 3-n: are the allele counts one column for each pool or population


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

Performs a simple least squares linear regression to predict missing allele counts per window:

Yₚ = XₚB

→ B̂ = inverse(XₚᵀXₚ) (XₚᵀYₚ).

Imputation is achieved by predicting the missing allele counts:

Ŷₘ = XₘB̂.

Where:

- Yₚ is the matrix of allele counts of pools with missing data at the loci without missing data (dimensions: mₚ non-missing loci × 7 alleles, nₘ pools with missing loci);
- Xₚ is the matrix of allele counts of pools without missing data at the loci without missing data in the other pools (dimensions: mₚ non-missing loci × 7 alleles, nₚ pools without missing loci);
- B̂ is the matrix of estimates of the effects of each pool without missing data on the allele counts of the pools with missing data (dimensions: nₚ pools without missing loci, nₘ pools with missing loci);
- inverse() is the Moore-Penrose pseudoinverse if the automatic Julia solver fails;
- Ŷₘ is the matrix of imputed allele counts of pools with missing data (dimensions: mₘ missing loci × 7 alleles, nₘ pools with missing loci); and
- Xₘ is the matrix of allele counts of pools without missing data at the loci with missing data in the other pools (dimensions: mₘ non-missing loci × 7 alleles, nₚ pools without missing loci).

The imputed allele counts are averaged across the windows sliding one locus at a time.


# Author
- Jeff Paril (jeffersonparil@gmail.com; https://orcid.org/0000-0002-5693-4123)
...

"""
function impute(str_filename_input; n_int_window_size=10, n_flt_maximum_fraction_of_pools_with_missing=0.5, n_flt_maximum_fraction_of_loci_with_missing=0.5, str_filename_output="output-imputed.syncx", n_int_thread_count=2)
    ###################################################################
    ### TEST
    # str_filename_input = "test.pileup"
    # n_int_window_size = 10
    # n_flt_maximum_fraction_of_pools_with_missing = 0.5
    # n_flt_maximum_fraction_of_loci_with_missing = 0.5
    # str_filename_output = "output-imputed.syncx"
    # n_int_thread_count = 2
    ###################################################################
    ### Opening remark
    println("")
    println("####################################################################")
    println("PoPoolImpute: imputation of population- and pool-level genotype data")
    println("(version 0.0.1; release 2021/12/30)")
    println("####################################################################")
    println("Input parameters:")
    @show str_filename_input
    @show n_int_window_size
    @show n_flt_maximum_fraction_of_pools_with_missing
    @show n_flt_maximum_fraction_of_loci_with_missing
    @show str_filename_output
    @show n_int_thread_count
    println("####################################################################")
    ### Define the full path to the input and output files
    if dirname(str_filename_input) == ""
        str_filename_input = string(pwd(), "/", str_filename_input)
    end
    if dirname(str_filename_output) == ""
        str_filename_output = string(pwd(), "/", str_filename_output)
    end
    ### Count the number of loci
    println("Counting the total number of loci in the input pileup file.")
    @show n_int_total_loci = countlines(str_filename_input)
    ### Define the size of each chunk to cut up the input file
    n_int_chunk_count = n_int_thread_count
    n_int_chuck_size = Int(ceil(n_int_total_loci / n_int_chunk_count))
    ### Rectify chunk size if it is less than twice the n_int_window_size
    if n_int_chuck_size < 2*n_int_window_size
        n_int_chunk_count = Int(floor(n_int_total_loci / (2*n_int_window_size)))
        n_int_chuck_size = Int(ceil(n_int_total_loci / n_int_chunk_count))
    end
    println(string("The input file was split into: ", n_int_chunk_count, " chunks of size: ", n_int_chuck_size, " loci each (at most)."))
    ### Cut up the input file
    open(str_filename_input) do FILE
        for i in 1:n_int_chunk_count
            j = 0
            file = open(string(str_filename_input, "-CHUNK_", i), "w")
            while (j < n_int_chuck_size) & (!eof(FILE))
                j += 1
                line = readline(FILE)
                write(file, string(line, '\n'))
            end
            close(file)
        end
    end
    println("Imputing.")
    @show n_int_thread_count
    if n_int_thread_count > 1
        ### Parallel for-loop
        println("Multi-threaded imputation.")
        @time _ = @sync @showprogress @distributed for i in 1:n_int_chunk_count
            str_filename_chunk_input = string(str_filename_input, "-CHUNK_", i)
            str_filename_chunk_output = string(str_filename_output, "-CHUNK_", i)
            fun_single_threaded_imputation(str_filename_chunk_input,
                                n_int_window_size=n_int_window_size,
                                n_flt_maximum_fraction_of_pools_with_missing=n_flt_maximum_fraction_of_pools_with_missing,
                                n_flt_maximum_fraction_of_loci_with_missing=n_flt_maximum_fraction_of_loci_with_missing,
                                str_filename_output=str_filename_chunk_output)
        end
    else
        println("Single-threaded imputation.")
        i = 1
        fun_single_threaded_imputation(str_filename_chunk_input,
                                n_int_window_size=n_int_window_size,
                                n_flt_maximum_fraction_of_pools_with_missing=n_flt_maximum_fraction_of_pools_with_missing,
                                n_flt_maximum_fraction_of_loci_with_missing=n_flt_maximum_fraction_of_loci_with_missing,
                                str_filename_output=str_filename_chunk_output)
    end
    ### Concatenate chunks
    ### NOTE: The issue here is that the imputed frequencies in the first and last windows were averaged from less imputation data points
    open(str_filename_output, "w") do FILE_OUT
        for i in 1:n_int_chunk_count
            file = open(string(str_filename_output, "-CHUNK_", i), "r")
            while !eof(file)
                line = readline(file)
                write(FILE_OUT, string(line, '\n'))
            end
            close(file)
            ### Clean-up
            rm(string(str_filename_input, "-CHUNK_", i))
            rm(string(str_filename_output, "-CHUNK_", i))
        end
    end
    ### Closing remark
    println("")
    println("####################################################################")
    println("Imputation successful. Please find the output file:")
    println(str_filename_output)
    println("####################################################################")
    ### Return code 0 for no error
    return(0)
end

end