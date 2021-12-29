module PoPoolImpute

using LinearAlgebra
include("functions.jl") ### load the functions with the module name qualifier
using .functions: fun_ascii_allele_states_to_counts_per_locus,
                                fun_ascii_allele_states_to_counts_per_window, 
                                fun_impute_per_window, 
                                fun_simple_progress_bar,
                                fun_writeout_inrun ### move the functions into scope, i.e. on need to use the module name qualifier

"""
# __________________________________________________________________
# PopoolImpute: impute allele frequency data from pool sequencing

`PopPoolImpute(str_filename_input; n_int_window_size=10, n_flt_maximum_fraction_of_pools_with_missing=0.5, n_flt_maximum_fraction_of_loci_with_missing=0.5, str_filename_output="output-imputed.syncx")`


# Inputs
1. *str_filename_input* [String]: filename of the genotype data in [pileup format (.pileup)](http://samtools.sourceforge.net/pileup.shtml)
1. *n_int_window_size* [Integer; default=10]: size of the sliding window across which the
...

# Outputs
1. *str_filename_output* [String; default="output-imputed.syncx"; comma-separated (.csv) file]
...

# Examples
```
str_filename_input = "test/test.pileup"
using PoPoolImpute
@time PoPoolImpute(str_filename_input)
```
# Details
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
    ### number of loci, alleles, and pools
    n_int_total_loci = countlines(str_filename_input)
    vec_allele_names=["A", "T", "C", "G", "INS", "DEL", "N"]
    n_int_allele_count = length(vec_allele_names)
    ### check if we have the expected number of columns in the first line of the pileup file, i.e. each pool has 3 columns each
    vec_line = split(readline(str_filename_input), "\t")
    n_int_pool_count = (length(vec_line) - 3) / 3
    if n_int_pool_count == round(n_int_pool_count)
        n_int_pool_count = Int(n_int_pool_count)
    else
        println("Ooopsss! Pileup file is corrupted.")
        println(string("Expected: ", Int(round(n_int_pool_count)), " pools but got ", n_int_pool_count, " pools instead."))
        println("Please check that each pool or population has 3 columns representing the depth, allele state, and allele quality.")
    end

    vec_str_NAME_OF_CHROMOSOME_OR_SCAFFOLD = []
    vec_int_POSITION = []

    n_int_counter_load_first_n_windows_lines_withe_readline = 0
    vec_str_input = []
    FILE = open(str_filename_input)

    for n_int_start_locus in 1:((n_int_total_loci-n_int_window_size) + 1)
        # n_int_start_locus = 25
        # @show n_int_start_locus
        fun_simple_progress_bar(n_int_start_locus, (n_int_total_loci-n_int_window_size) + 1, 50)
        # vec_str_input = readlines(str_filename_input)[n_int_start_locus:(n_int_start_locus+n_int_window_size-1)]
        ### if we already have n_int_window_size lines then just remove the old locus and replace with the next one since we have sliding windows
        if n_int_counter_load_first_n_windows_lines_withe_readline == n_int_window_size
            vec_str_input[1:(end-1)] = vec_str_input[2:end]
            vec_str_input[end] = readline(FILE)
        else
            ### parse n_int_window_size lines
            while n_int_counter_load_first_n_windows_lines_withe_readline < n_int_window_size
                push!(vec_str_input, readline(FILE))
                n_int_counter_load_first_n_windows_lines_withe_readline += 1
            end
        end
        ### parse pileup into a matrix of allele + locus IDs ### current bottleneck O(n) I think - will have to properly measure...
        vec_str_name_of_chromosome_or_scaffold, vec_int_position, mat_int_window_counts = fun_ascii_allele_states_to_counts_per_window(vec_str_input, vec_allele_names)
        n_bool_window_with_at_least_one_missing_locus = sum(ismissing.(mat_int_window_counts)) > 0
        ### If we have missing loci then impute, else just add the allele counts without missing information
        if n_bool_window_with_at_least_one_missing_locus
            ### impute
            mat_imputed, vec_bool_idx_pools_with_missing_loci, vec_bool_idx_loci_missing = fun_impute_per_window(mat_int_window_counts, n_flt_maximum_fraction_of_pools_with_missing, n_flt_maximum_fraction_of_loci_with_missing)
            ### replace missing data with the imputed or predicted allele counts if we were able to impute, i.e. we got at mot most 50% of the pools with missing loci, and 50% of the loci with missing data
            if !ismissing(mat_imputed)
                mat_int_window_counts[vec_bool_idx_loci_missing, vec_bool_idx_pools_with_missing_loci] = mat_imputed
            end
        else
            mat_imputed = "Not missing but no imputation needed since no loci were missing."
        end
        ### add allele counts imputed and non-missing + loci information
        if !ismissing(mat_imputed)
            if !@isdefined mat_int_ALLELE_COUNTS
                ### initialise the output matrix of allele counts and append the chromosom or scaffold names and the position
                global mat_int_ALLELE_COUNTS = Int.(mat_int_window_counts)
                append!(vec_str_NAME_OF_CHROMOSOME_OR_SCAFFOLD, vec_str_name_of_chromosome_or_scaffold)
                append!(vec_int_POSITION, vec_int_position)
            else
                ### append loci information into the output matrix and vectors
                vec_bool_idx_loci_existing_loci = [x ∈ vec_str_NAME_OF_CHROMOSOME_OR_SCAFFOLD for x in vec_str_name_of_chromosome_or_scaffold] .&
                                                  [x ∈ vec_int_POSITION for x in vec_int_position]
                vec_bool_idx_loci_new_loci_to_add = .!vec_bool_idx_loci_existing_loci
                # @show vec_int_POSITION
                # @show vec_int_position
                
                ### Do we have more than 1 new loci (happens if we skip loci due to the inability to impute because of too many missing data)
                n_int_new_loci_count = sum(vec_bool_idx_loci_new_loci_to_add)
                # @show n_int_new_loci_count

                ### if we have imputed allele counts, then use the average of the imputed allele counts
                if mat_imputed != "Not missing but no imputation needed since no loci were missing."
                    if n_int_new_loci_count < n_int_window_size
                        ### Compute the average allele count
                        vec_idx_bool_loci_missing_less_new_locus = vec_bool_idx_loci_missing[1:(end-(n_int_allele_count*n_int_new_loci_count))]
                        mat_int_allele_counts_tail_end_old = mat_int_ALLELE_COUNTS[(end-(n_int_allele_count*(n_int_window_size-n_int_new_loci_count))+1):end, :]
                        mat_int_allele_counts_tail_end_new = mat_int_window_counts[1:(end-(n_int_allele_count*n_int_new_loci_count)), :]
                        mat_bool_idx_loci_missing_less_new_locus = reshape(vec_idx_bool_loci_missing_less_new_locus, (n_int_allele_count, n_int_window_size-n_int_new_loci_count))'
                        vec_int_imputed_loci_counter = ones(Int, n_int_window_size-n_int_new_loci_count)
                        vec_bool_index_for_vec_int_imputed_loci_counter = (sum(mat_bool_idx_loci_missing_less_new_locus, dims=2) .> 0)[:,1]
                        vec_n = repeat(vec_int_imputed_loci_counter[vec_bool_index_for_vec_int_imputed_loci_counter], inner=n_int_allele_count)
                        A = mat_int_allele_counts_tail_end_old[vec_idx_bool_loci_missing_less_new_locus, :]
                        B = mat_int_allele_counts_tail_end_new[vec_idx_bool_loci_missing_less_new_locus, :]
                        C = Int.(round.( ( (A.+(B./vec_n)) ./ 2 ) .* ( (2 .* vec_n) ./ (vec_n .+ 1) ) ))
                        ### Update allele counts with the average
                        mat_int_ALLELE_COUNTS[(end-(n_int_allele_count*(n_int_window_size-n_int_new_loci_count))+1):end, :][vec_idx_bool_loci_missing_less_new_locus, :] = C
                        ### Update counter
                        vec_int_imputed_loci_counter = vec_int_imputed_loci_counter .+ vec_bool_index_for_vec_int_imputed_loci_counter
                        vec_int_imputed_loci_counter[1:(end-1)] = vec_int_imputed_loci_counter[2:(end-0)]
                        vec_int_imputed_loci_counter[end] = 1
                    end
                end

                ### Save the file per window because it's nice to have the output written into disk rather than memory in case anything unsavory happens prior to finishing the entire job - then at least we'll have a partial output rather than nothing at all
                ### Save a locus once we're done with the trailing end of the previous window
                if (n_int_start_locus >= 2)
                    println("Debugging.............................. --- 0")
                    fun_writeout_inrun(vec_str_NAME_OF_CHROMOSOME_OR_SCAFFOLD[1],
                                       vec_int_POSITION[1],
                                       mat_int_ALLELE_COUNTS[1:n_int_allele_count, :],
                                       n_int_allele_count,
                                       str_filename_output)
                end
                
                ### Update with new loci keeping the size constant by removing the loci out of the window and adding the new loci entering the window
                if n_int_new_loci_count < n_int_window_size
                    mat_int_ALLELE_COUNTS[1:(end-(n_int_allele_count*n_int_new_loci_count)), :] = mat_int_ALLELE_COUNTS[((n_int_allele_count*n_int_new_loci_count)+1):end, :]
                    mat_int_ALLELE_COUNTS[((end-(n_int_allele_count*n_int_new_loci_count))+1):end, :] = mat_int_window_counts[repeat(vec_bool_idx_loci_new_loci_to_add, inner=n_int_allele_count), :]
                    vec_str_NAME_OF_CHROMOSOME_OR_SCAFFOLD[1:(end-n_int_new_loci_count)] = vec_str_NAME_OF_CHROMOSOME_OR_SCAFFOLD[(n_int_new_loci_count+1):end]
                    vec_str_NAME_OF_CHROMOSOME_OR_SCAFFOLD[((end-n_int_new_loci_count)+1):end] = vec_str_name_of_chromosome_or_scaffold[vec_bool_idx_loci_new_loci_to_add]
                    vec_int_POSITION[1:(end-n_int_new_loci_count)] = vec_int_POSITION[(n_int_new_loci_count+1):end]
                    vec_int_POSITION[((end-n_int_new_loci_count)+1):end] = vec_int_position[vec_bool_idx_loci_new_loci_to_add]
                else
                    mat_int_ALLELE_COUNTS[1:end,:] = mat_int_window_counts
                    vec_str_NAME_OF_CHROMOSOME_OR_SCAFFOLD[1:end] = vec_str_name_of_chromosome_or_scaffold
                    vec_int_POSITION[1:end] = vec_int_position
                end
                ### If we reach the end of the file (1 window offset)
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
    ### Return zero is all is well.
    return(0)
end

end