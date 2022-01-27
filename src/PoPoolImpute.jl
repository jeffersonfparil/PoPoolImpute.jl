module PoPoolImpute

using ProgressMeter
using Distributed
using Dates ### included in base julia installation as part of the standard library and so no need to install as a dependency in .github/workflows/julia.yml
### Load the functions, and move them into scope
include("functions.jl")
using .functions: fun_ascii_allele_states_to_counts_per_locus,
                                fun_ascii_allele_states_to_counts_per_window,
                                fun_impute_per_window,
                                fun_simple_progress_bar,
                                fun_split_pileup,
                                fun_writeout_inrun,
                                fun_single_threaded_imputation
### Documentation
"""
# ____________________________________________________________________
# PoPoolImpute: imputation of population- and pool-level genotype data


# Usage
`PopPoolImpute.impute(str_filename_input; n_int_window_size=10, n_flt_maximum_fraction_of_pools_with_missing=0.5, n_flt_maximum_fraction_of_loci_with_missing=0.5, str_filename_output="output-imputed.syncx", bool_OLS_dist=false, n_int_thread_count=1)`


# Inputs
1. str\\_filename\\_input [String]: filename of the genotype data in [pileup format (.pileup)](http://samtools.sourceforge.net/pileup.shtml)
2. n\\_int\\_window\\_size [Integer; default=10]: size of the sliding window across which imputation is performed
3. n\\_flt\\_maximum\\_fraction\\_of\\_pools\\_with\\_missing [Float; default=0.5]: maximum tolerable fraction of the pools with at least one missing locus
4. n\\_flt\\_maximum\\_fraction\\_of\\_loci\\_with\\_missing [Float; default=0.5]: maximum tolerable fraction of the loci with missing data
5. bool\\_OLS\\_dist [Boolean; default=false]: perform ordinary least squares regression with pairwise loci distances covariate
6. n\\_int\\_thread\\_count [Integer; default=1]: number of computing threads to use


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
function impute(str_filename_input; n_int_window_size=10, n_flt_maximum_fraction_of_pools_with_missing=0.5, n_flt_maximum_fraction_of_loci_with_missing=0.5, str_filename_output="output-imputed.syncx", bool_OLS_dist=false, n_int_thread_count=2)
    ###################################################################
    ### TEST
    # cd("/home/jeff/Documents/PoPoolImpute.jl/test")
    # str_filename_input = "out_simissing.pileup"
    # n_int_window_size = 10
    # n_flt_maximum_fraction_of_pools_with_missing = 0.5
    # n_flt_maximum_fraction_of_loci_with_missing = 0.5
    # str_filename_output = "output-imputed.syncx"
    # n_int_thread_count = 3
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
    println(string("Time: ", Dates.format(now(), "Y-u-dd(E)THH:MM")))
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
    ### Split the input pileup file if we can afford parallel processing, i.e. (n_int_thread_count > 1) since (n_int_thread_count = n_int_chunk_count)
    if n_int_chunk_count > 1
        fun_split_pileup(str_filename_input,
                         n_int_chunk_count=n_int_chunk_count,
                         n_int_chuck_size=n_int_chuck_size,
                         n_int_window_size=n_int_window_size,
                         n_bool_add_leading_trailing_windows=true)
    end
    println("Imputing.")
    @show n_int_thread_count
    if n_int_thread_count > 1
        ### Parallel for-loop for each chunk of the split input file
        println("Multi-threaded imputation.")
        @time _ = @sync @showprogress @distributed for i in 1:n_int_chunk_count
            str_filename_chunk_input = string(str_filename_input, "-CHUNK_", i)
            str_filename_chunk_output = string(str_filename_output, "-CHUNK_", i)
            ### Trim-out the leading and/or trailing windows
            n_bool_skip_leading_window = i>1 ### first chunk does not have a leading window
            n_bool_skip_trailing_window = i<n_int_chunk_count ### last chunk does not have trailing window
            fun_single_threaded_imputation(str_filename_chunk_input,
                                n_int_window_size=n_int_window_size,
                                n_flt_maximum_fraction_of_pools_with_missing=n_flt_maximum_fraction_of_pools_with_missing,
                                n_flt_maximum_fraction_of_loci_with_missing=n_flt_maximum_fraction_of_loci_with_missing,
                                str_filename_output=str_filename_chunk_output,
                                n_bool_skip_leading_window=n_bool_skip_leading_window,
                                n_bool_skip_trailing_window=n_bool_skip_trailing_window,
                                bool_OLS_dist=bool_OLS_dist)
        end
    else
        ### Non-parallel imputation for unsplit input file
        println("Single-threaded imputation.")
        fun_single_threaded_imputation(str_filename_input,
                                n_int_window_size=n_int_window_size,
                                n_flt_maximum_fraction_of_pools_with_missing=n_flt_maximum_fraction_of_pools_with_missing,
                                n_flt_maximum_fraction_of_loci_with_missing=n_flt_maximum_fraction_of_loci_with_missing,
                                str_filename_output=str_filename_output,
                                n_bool_skip_leading_window=false,
                                n_bool_skip_trailing_window=false,
                                bool_OLS_dist=bool_OLS_dist)
    end
    ### Concatenate chunks
    if n_int_chunk_count > 1
        open(str_filename_output, "w") do FILE_OUT
            for i in 1:n_int_chunk_count
                str_filename_chunk = string(str_filename_output, "-CHUNK_", i)
                ### Check if the chunk exists (the chunk file may not exist if no loci were kept)
                if isfile(str_filename_chunk)
                    file = open(str_filename_chunk, "r")
                    while !eof(file)
                        write(FILE_OUT, string(readline(file), '\n'))
                    end
                    close(file)
                    ### Clean-up
                    rm(string(str_filename_output, "-CHUNK_", i))
                else
                    ### Message when a chunk was not imputed at all because the data is too sparse,
                    ### where sparsity is discatated by the maximum number of loci and pools allowed to be missing
                    ### as well as successful parameter estimation of the linear regression model.
                    println("No imputation output for: ")
                    println(string(str_filename_input, "-CHUNK_", i))
                end
                ### Clean-up
                rm(string(str_filename_input, "-CHUNK_", i))
            end
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