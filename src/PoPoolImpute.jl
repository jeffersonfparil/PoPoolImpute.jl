module PoPoolImpute

### Load linear algebra library for the Moore-Penrose pseudoinverse if the automatic solver fails
using LinearAlgebra
### Load the functions, and move them into scope
include("functions.jl")
using .functions: fun_ascii_allele_states_to_counts_per_locus,
                                fun_ascii_allele_states_to_counts_per_window, 
                                fun_impute_per_window, 
                                fun_simple_progress_bar,
                                fun_writeout_inrun
### Documentation
"""
# ____________________________________________________________________
# PoPoolImpute: imputation of population- and pool-level genotype data

`PopPoolImpute(str_filename_input; n_int_window_size=10, n_flt_maximum_fraction_of_pools_with_missing=0.5, n_flt_maximum_fraction_of_loci_with_missing=0.5, str_filename_output="output-imputed.syncx")`


# Inputs
1. str\\_filename\\_input [String]: filename of the genotype data in [pileup format (.pileup)](http://samtools.sourceforge.net/pileup.shtml)
2. n\\_int\\_window\\_size [Integer; default=10]: size of the sliding window across which imputation is performed
2. n\\_flt\\_maximum\\_fraction\\_of\\_pools\\_with\\_missing [Float; default=0.5]: maximum tolerable fraction of the pools with at least one missing locus
2. n\\_flt\\_maximum\\_fraction\\_of\\_loci\\_with\\_missing [Float; default=0.5]: maximum tolerable fraction of the loci with missing data
...

# Output
str\\_filename\\_output [String; default="output-imputed.syncx"; comma-separated file]

Syncx format (after popoolation2's sync or synchronised pileup file format):
- Column 1:   chromosome or scaffold name
- Column 2:   locus position repeated 7 times corresponding to alleles "A", "T", "C", "G", "INS", "DEL", "N", where "INS" is insertion, "DEL" is deletion, and "N" is unclassified
- Column 3-n: are the allele counts one column for each pool or population
...

# Examples
```
using PoPoolImpute
str_filename_input = "test.pileup"
PoPoolImpute(str_filename_input)
PoPoolImpute(str_filename_input, n_int_window_size=20, str_filename_output="test-2.syncx")
PoPoolImpute(str_filename_input, n_flt_maximum_fraction_of_pools_with_missing=0.2, str_filename_output="test-3.syncx")
```
# Details

Performs a simple least squares linear regression to predict missing allele counts per window:

Yₚ = XₚB

→ B̂ = inverse(Xₚ'Xₚ) (Xₚ'Yₚ).

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

# Authors
- Jeff Paril (jeffersonparil@gmail.com; https://orcid.org/0000-0002-5693-4123)
...

"""
function PopPoolImpute(str_filename_input; n_int_window_size=10, n_flt_maximum_fraction_of_pools_with_missing=0.5, n_flt_maximum_fraction_of_loci_with_missing=0.5, str_filename_output="output-imputed.syncx")
    ###################################################################
    ### TEST
    # cd("/home/jeff/Documents/PoPoolImpute.jl/test")
    # str_filename_input = "out_simissing.pileup"
    # n_int_window_size = 10
    # n_flt_maximum_fraction_of_pools_with_missing = 0.5
    # n_flt_maximum_fraction_of_loci_with_missing = 0.5
    # str_filename_output = "output-imputed.syncx"
    ###################################################################
    ### Opening remark
    println("")
    println("####################################################################")
    println("PoPoolImpute: imputation of population- and pool-level genotype data")
    println("(version 0.0.1; release 2021/12/30)")
    println("####################################################################")
    ### Count the number of loci, alleles, and pools (check if we have the expected number of columns in the first line of the pileup file, i.e. each pool has 3 columns each)
    n_int_total_loci = countlines(str_filename_input)
    vec_allele_names=["A", "T", "C", "G", "INS", "DEL", "N"]
    n_int_allele_count = length(vec_allele_names)
    FILE_to_find_pool_count = open(str_filename_input)
    vec_line = split(readline(FILE_to_find_pool_count), "\t")
    close(FILE_to_find_pool_count)
    n_int_pool_count = (length(vec_line) - 3) / 3
    if n_int_pool_count == round(n_int_pool_count)
        n_int_pool_count = Int(n_int_pool_count)
    else
        println("Ooopsss! Pileup file is corrupted.")
        println(string("Expected: ", Int(round(n_int_pool_count)), " pools but got ", n_int_pool_count, " pools instead."))
        println("Please check that each pool or population has 3 columns representing the depth, allele state, and allele quality.")
    end
    ### Initialise the vectors of scaffold or chromosome names and positions
    vec_str_NAME_OF_CHROMOSOME_OR_SCAFFOLD = []
    vec_int_POSITION = []
    mat_int_ALLELE_COUNTS = nothing
    ### Initialise loci within window counter, vector containing each line (will fit a maximum of n_int_window_size loci), and open the file
    n_int_counter_load_first_n_windows_lines_withe_readline = 0
    vec_str_input = []
    FILE = open(str_filename_input)
    ### Iterate per line until we reach the first loci of the last sliding window
    for n_int_start_locus in 1:((n_int_total_loci-n_int_window_size) + 1)
        # n_int_start_locus = 25
        # @show n_int_start_locus
        fun_simple_progress_bar(n_int_start_locus, (n_int_total_loci-n_int_window_size) + 1, 50)
        ### If we already have "n_int_window_size" lines then just remove the old locus and replace with the next one since we have sliding windows (sliding one locus at a time)
        if n_int_counter_load_first_n_windows_lines_withe_readline == n_int_window_size
            vec_str_input[1:(end-1)] = vec_str_input[2:end]
            vec_str_input[end] = readline(FILE) ### read the next line of the file (iterates per line as long as "close(FILE)" has not been executed)
        else
            ### Fill "vec_str_input" so we have "n_int_window_size" loci
            while n_int_counter_load_first_n_windows_lines_withe_readline < n_int_window_size
                push!(vec_str_input, readline(FILE))
                n_int_counter_load_first_n_windows_lines_withe_readline += 1
            end
        end
        ### Parse each window (contained in "vec_str_input") to extract the chromosome or scaffold names, positions, and the matrix of allele counts with "(n_int_window_size*vec_allele_names) x n_int_pool_count" dimensions
        vec_str_name_of_chromosome_or_scaffold, vec_int_position, mat_int_window_counts = fun_ascii_allele_states_to_counts_per_window(vec_str_input, vec_allele_names)
        n_bool_window_with_at_least_one_missing_locus = sum(ismissing.(mat_int_window_counts)) > 0
        ### If we have missing loci then impute, else just add the allele counts without missing information
        if n_bool_window_with_at_least_one_missing_locus
            ### Impute by regressing allele counts of the pools with missing data against the allele counts of the pools without missing data in the window; and then predict the missing allele counts
            mat_imputed, vec_bool_idx_pools_with_missing_loci, vec_bool_idx_loci_missing = fun_impute_per_window(mat_int_window_counts, n_flt_maximum_fraction_of_pools_with_missing, n_flt_maximum_fraction_of_loci_with_missing)
            ### Replace missing data with the imputed allele counts if we were able to impute, i.e. we got at mot most "n_flt_maximum_fraction_of_pools_with_missing" of the pools with missing loci, and "n_flt_maximum_fraction_of_loci_with_missing" of the loci with missing data
            if !ismissing(mat_imputed)
                mat_int_window_counts[vec_bool_idx_loci_missing, vec_bool_idx_pools_with_missing_loci] = mat_imputed
            end
        else
            mat_imputed = "Not missing but no imputation needed since no loci were missing."
        end
        ### Add allele counts imputed and non-missing + loci information
        if !ismissing(mat_imputed)
            if isnothing(mat_int_ALLELE_COUNTS)
                ### initialise the output matrix of allele counts and append the chromosom or scaffold names and the position
                mat_int_ALLELE_COUNTS = Int.(mat_int_window_counts)
                append!(vec_str_NAME_OF_CHROMOSOME_OR_SCAFFOLD, vec_str_name_of_chromosome_or_scaffold)
                append!(vec_int_POSITION, vec_int_position)
            else
                ### append loci information into the output matrix and vectors
                vec_bool_idx_loci_existing_loci = [x ∈ vec_str_NAME_OF_CHROMOSOME_OR_SCAFFOLD for x in vec_str_name_of_chromosome_or_scaffold] .&
                                                  [x ∈ vec_int_POSITION for x in vec_int_position]
                vec_bool_idx_loci_new_loci_to_add = .!vec_bool_idx_loci_existing_loci
                ### Do we have more than 1 new loci (happens if we skip loci due to the inability to impute because of too many missing data)
                n_int_new_loci_count = sum(vec_bool_idx_loci_new_loci_to_add)
                ### if we have imputed allele counts, then use the average of the imputed allele counts
                if mat_imputed != "Not missing but no imputation needed since no loci were missing."
                    if n_int_new_loci_count < n_int_window_size
                        ### Compute the average imputed allele counts by updating the average given new imputed allele counts
                        vec_idx_bool_loci_missing_less_new_loci = vec_bool_idx_loci_missing[1:(end-(n_int_allele_count*n_int_new_loci_count))]
                        mat_int_allele_counts_tail_end_old = mat_int_ALLELE_COUNTS[(end-(n_int_allele_count*(n_int_window_size-n_int_new_loci_count))+1):end, :]
                        mat_int_allele_counts_tail_end_new = mat_int_window_counts[1:(end-(n_int_allele_count*n_int_new_loci_count)), :]
                        mat_bool_idx_loci_missing_less_new_locus = reshape(vec_idx_bool_loci_missing_less_new_loci, (n_int_allele_count, n_int_window_size-n_int_new_loci_count))'
                        vec_int_imputed_loci_counter = ones(Int, n_int_window_size-n_int_new_loci_count)
                        vec_bool_index_for_vec_int_imputed_loci_counter = (sum(mat_bool_idx_loci_missing_less_new_locus, dims=2) .> 0)[:,1]
                        vec_n = repeat(vec_int_imputed_loci_counter[vec_bool_index_for_vec_int_imputed_loci_counter], inner=n_int_allele_count)
                        A = mat_int_allele_counts_tail_end_old[vec_idx_bool_loci_missing_less_new_loci, :]
                        B = mat_int_allele_counts_tail_end_new[vec_idx_bool_loci_missing_less_new_loci, :]
                        C = Int.(round.( ( (A.+(B./vec_n)) ./ 2 ) .* ( (2 .* vec_n) ./ (vec_n .+ 1) ) ))
                        ### Update allele counts with the average
                        mat_int_ALLELE_COUNTS[(end-(n_int_allele_count*(n_int_window_size-n_int_new_loci_count))+1):end, :][vec_idx_bool_loci_missing_less_new_loci, :] = C
                        ### Update the counter which we use to compute the updated average allele count
                        vec_int_imputed_loci_counter = vec_int_imputed_loci_counter .+ vec_bool_index_for_vec_int_imputed_loci_counter
                        vec_int_imputed_loci_counter[1:(end-1)] = vec_int_imputed_loci_counter[2:(end-0)]
                        vec_int_imputed_loci_counter[end] = 1
                    end
                end
                ### Save the file per window because it's nice to have the output written into disk rather than memory in case anything unsavory happens prior to finishing the entire job - then at least we'll have a partial output rather than nothing at all
                if (n_int_start_locus >= 2)
                    ### Save a locus once we're done with the trailing end of the previous window
                    fun_writeout_inrun(vec_str_NAME_OF_CHROMOSOME_OR_SCAFFOLD[1],
                                       vec_int_POSITION[1],
                                       mat_int_ALLELE_COUNTS[1:n_int_allele_count, :],
                                       n_int_allele_count,
                                       str_filename_output)
                end
                ### Update the matrix of allele counts, and vectors of loci coordinates with the new loci keeping the size constant by removing the loci out of the window and adding the new loci entering the window
                if n_int_new_loci_count < n_int_window_size
                    mat_int_ALLELE_COUNTS[1:(end-(n_int_allele_count*n_int_new_loci_count)), :] = mat_int_ALLELE_COUNTS[((n_int_allele_count*n_int_new_loci_count)+1):end, :]
                    mat_int_ALLELE_COUNTS[((end-(n_int_allele_count*n_int_new_loci_count))+1):end, :] = mat_int_window_counts[repeat(vec_bool_idx_loci_new_loci_to_add, inner=n_int_allele_count), :]
                    vec_str_NAME_OF_CHROMOSOME_OR_SCAFFOLD[1:(end-n_int_new_loci_count)] = vec_str_NAME_OF_CHROMOSOME_OR_SCAFFOLD[(n_int_new_loci_count+1):end]
                    vec_str_NAME_OF_CHROMOSOME_OR_SCAFFOLD[((end-n_int_new_loci_count)+1):end] = vec_str_name_of_chromosome_or_scaffold[vec_bool_idx_loci_new_loci_to_add]
                    vec_int_POSITION[1:(end-n_int_new_loci_count)] = vec_int_POSITION[(n_int_new_loci_count+1):end]
                    vec_int_POSITION[((end-n_int_new_loci_count)+1):end] = vec_int_position[vec_bool_idx_loci_new_loci_to_add]
                else
                    ### If none of the old and new windows overlap, i.e. when "mat_int_window_counts" loci were skipped because of too much missing data which rendered imputation impossible,
                    ### then replace the matrix of allele counts, and vectors of loci coordinates with the new loci
                    mat_int_ALLELE_COUNTS[1:end,:] = mat_int_window_counts
                    vec_str_NAME_OF_CHROMOSOME_OR_SCAFFOLD[1:end] = vec_str_name_of_chromosome_or_scaffold
                    vec_int_POSITION[1:end] = vec_int_position
                end
                ### If we reach the end of the file offset by one window then save the last window
                if n_int_start_locus == ((n_int_total_loci-n_int_window_size) + 1)
                    for i in 1:length(vec_str_NAME_OF_CHROMOSOME_OR_SCAFFOLD)
                        fun_writeout_inrun(vec_str_NAME_OF_CHROMOSOME_OR_SCAFFOLD[i],
                                        vec_int_POSITION[i],
                                        mat_int_ALLELE_COUNTS[(((i-1)*n_int_allele_count)+1):(i*n_int_allele_count), :],
                                        n_int_allele_count,
                                        str_filename_output)
                    end
                end
            end
        end
    end
    ### Note: loci which cannot be imputed are not included in the output file
    ### Close input file
    close(FILE)
    ### Message to indicate the output file
    println("")
    println("####################################################################")
    println("Imputation successful. Please find the output file:")
    println(str_filename_output)
    println("####################################################################")
    ### Return zero to indicate that all is well
    return(0)
end

end