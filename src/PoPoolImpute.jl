module PoPoolImpute

include("functions.jl") ### load the functions with the module name qualifier
using .functions: fun_ascii_allele_states_to_counts_per_locus,
                                fun_ascii_allele_states_to_counts_per_window, 
                                fun_impute_per_window, 
                                fun_simple_progress_bar,
                                fun_writeout_inrun ### move the functions into scope, i.e. on need to use the module name qualifier

"""
# __________________________________________________________________
# PopoolImpute: impute allele frequency data from pool sequencing

`PopPoolImpute(str_filename_input; vec_allele_names=["A", "T", "C", "G", "INS", "DEL", "N"], n_int_window_size=10, str_filename_output="output-imputed.syncx")`


# Inputs
1. *str_filename_input* [String]: filename of the genotype data in [pileup format (.sync)](http://samtools.sourceforge.net/pileup.shtml)
...

# Outputs
1. *str_filename_output* [String; default="output-imputed.syncx"; comma-separated (.csv) file]
...

# Examples
```
### Input files:
str_filename_input = "test/test.pileup"
### Single-threaded execution:
using PoPoolImpute
@time PoPoolImpute(str_filename_input)
### Multi-threaded execution (parallel execution only applicable to model=="GWAlpha"):
using Distributed
Distributed.addprocs(length(Sys.cpu_info())-1)
@everywhere using PoPoolImpute
@time PoPoolImpute(str_filename_input)
```
# Details
...

"""
function PopPoolImpute(str_filename_input; vec_allele_names=["A", "T", "C", "G", "INS", "DEL", "N"], n_int_window_size=10, str_filename_output="output-imputed.syncx")
    ###################################################################
    ### TEST
    # str_filename_input = "out_simissing.mpileup.debugging"
    # vec_allele_names = ["A", "T", "C", "G", "INS", "DEL", "N"]
    # n_int_window_size = 10
    # str_filename_output = "output-imputed.syncx"
    ###################################################################
    ### number of loci, alleles, and pools
    n_int_total_loci = countlines(str_filename_input)
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
    vec_int_imputed_loci_counter = ones(Int, n_int_window_size-1)

    n_int_counter_load_first_n_windows_lines_withe_readline = 0
    vec_str_input = []
    FILE = open(str_filename_input)

    for n_int_start_locus in 1:((n_int_total_loci-n_int_window_size) + 1)
        # n_int_start_locus = 2
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
            mat_imputed, vec_bool_idx_pools_with_missing_loci, vec_bool_idx_loci_missing = fun_impute_per_window(mat_int_window_counts)
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
                vec_bool_idx_loci_existing_loci = [x ∈ vec_str_NAME_OF_CHROMOSOME_OR_SCAFFOLD[(end-(n_int_window_size-1)):end] for x in vec_str_name_of_chromosome_or_scaffold] .&
                                                  [x ∈ vec_int_POSITION[(end-(n_int_window_size-1)):end] for x in vec_int_position]
                vec_bool_idx_loci_new_loci_to_add = .!vec_bool_idx_loci_existing_loci
                ### if we have imputed allele counts, then use the average of the imputed allele counts
                if mat_imputed != "Not missing but no imputation needed since no loci were missing."
                    ### Compute the average allele count
                    vec_idx_bool_loci_missing_less_new_locus = vec_bool_idx_loci_missing[1:(end-n_int_allele_count)]
                    mat_int_allele_counts_tail_end_old = mat_int_ALLELE_COUNTS[(end-(n_int_allele_count*(n_int_window_size-1))+1):end, :]
                    mat_int_allele_counts_tail_end_new = mat_int_window_counts[1:(end-n_int_allele_count), :]
                    mat_bool_idx_loci_missing_less_new_locus = reshape(vec_idx_bool_loci_missing_less_new_locus, (n_int_allele_count, n_int_window_size-1))'
                    vec_bool_index_for_vec_int_imputed_loci_counter = (sum(mat_bool_idx_loci_missing_less_new_locus, dims=2) .> 0)[:,1]
                    vec_n = repeat(vec_int_imputed_loci_counter[vec_bool_index_for_vec_int_imputed_loci_counter], inner=n_int_allele_count)
                    # @show n_int_start_locus
                    # @show vec_n
                    # println("Debugging..............................")
                    A = mat_int_allele_counts_tail_end_old[vec_idx_bool_loci_missing_less_new_locus, :]
                    B = mat_int_allele_counts_tail_end_new[vec_idx_bool_loci_missing_less_new_locus, :]
                    C = Int.(round.( ( (A.+(B./vec_n)) ./ 2 ) .* ( (2 .* vec_n) ./ (vec_n .+ 1) ) ))
                    ### Update allele counts with the average
                    mat_int_ALLELE_COUNTS[(end-(n_int_allele_count*(n_int_window_size-1))+1):end, :][vec_idx_bool_loci_missing_less_new_locus, :] = C
                    ### Update counter
                    vec_int_imputed_loci_counter = vec_int_imputed_loci_counter .+ vec_bool_index_for_vec_int_imputed_loci_counter
                    vec_int_imputed_loci_counter[1:(end-1)] = vec_int_imputed_loci_counter[2:(end-0)]
                    vec_int_imputed_loci_counter[end] = 1
                end
                ### add new loci regardless of whether or not we have imputed missing allele counts
                # mat_int_ALLELE_COUNTS = vcat(mat_int_ALLELE_COUNTS, mat_int_window_counts[repeat(vec_bool_idx_loci_new_loci_to_add, inner=n_int_allele_count), :])
                # append!(vec_str_NAME_OF_CHROMOSOME_OR_SCAFFOLD, vec_str_name_of_chromosome_or_scaffold[vec_bool_idx_loci_new_loci_to_add])
                # append!(vec_int_POSITION, vec_int_position[vec_bool_idx_loci_new_loci_to_add])

                ### Save the file per window because it's nice to have the output written into disk rather than memory in case anything unsavory happens prior to finishing the entire job - then at least we'll have a partial output rather than nothing at all
                ### Save a locus once we're done with the trailing end of the previous window
                if (n_int_start_locus >= 2)
                    # OUT = hcat(repeat([vec_str_NAME_OF_CHROMOSOME_OR_SCAFFOLD[end-n_int_window_size]], inner=n_int_allele_count),
                    #            repeat([vec_int_POSITION[end-n_int_window_size]], inner=n_int_allele_count),
                    #            mat_int_ALLELE_COUNTS[((end-(n_int_window_size*n_int_allele_count))-n_int_allele_count+1):(end-(n_int_window_size*n_int_allele_count)), :]
                    #           )
                    # OUT = hcat(repeat([vec_str_NAME_OF_CHROMOSOME_OR_SCAFFOLD[1]], inner=n_int_allele_count),
                    #            repeat([vec_int_POSITION[1]], inner=n_int_allele_count),
                    #            mat_int_ALLELE_COUNTS[((end-(n_int_window_size*n_int_allele_count))+1):((end-(n_int_window_size*n_int_allele_count))+n_int_allele_count), :]
                    #           )
                    # out = join([join(x,',') for x in eachrow(OUT)], '\n')
                    # file = open(str_filename_output, "a")
                    # write(file, string(out, '\n'))
                    # close(file)
                    fun_writeout_inrun(vec_str_NAME_OF_CHROMOSOME_OR_SCAFFOLD[1],
                                    vec_int_POSITION[1],
                                    mat_int_ALLELE_COUNTS[((end-(n_int_window_size*n_int_allele_count))+1):((end-(n_int_window_size*n_int_allele_count))+n_int_allele_count), :],
                                    n_int_allele_count,
                                    str_filename_output)
                end
                
                ### Update with new loci keeping the size constant by removing the loci out of the window and adding the new loci entering the window
                mat_int_ALLELE_COUNTS[1:(end-n_int_allele_count), :] = mat_int_ALLELE_COUNTS[(n_int_allele_count+1):end, :]
                mat_int_ALLELE_COUNTS[((end-n_int_allele_count)+1):end, :] = mat_int_window_counts[repeat(vec_bool_idx_loci_new_loci_to_add, inner=n_int_allele_count), :]
                vec_str_NAME_OF_CHROMOSOME_OR_SCAFFOLD[1:(end-1)] = vec_str_NAME_OF_CHROMOSOME_OR_SCAFFOLD[2:end]
                vec_str_NAME_OF_CHROMOSOME_OR_SCAFFOLD[end] = vec_str_name_of_chromosome_or_scaffold[vec_bool_idx_loci_new_loci_to_add][1]
                vec_int_POSITION[1:(end-1)] = vec_int_POSITION[2:end]
                vec_int_POSITION[end] = vec_int_position[vec_bool_idx_loci_new_loci_to_add][1]
                ### If we reach the end of the file (1 window offset)
                if n_int_start_locus == ((n_int_total_loci-n_int_window_size) + 1)
                    # OUT = hcat(repeat(vec_str_NAME_OF_CHROMOSOME_OR_SCAFFOLD[(end-n_int_window_size):end], inner=n_int_allele_count),
                    #            repeat(vec_int_POSITION[(end-n_int_window_size):end], inner=n_int_allele_count),
                    #            mat_int_ALLELE_COUNTS[((end-(n_int_window_size*n_int_allele_count))-n_int_allele_count+1):end, :]
                    #           )
                    # OUT = hcat(repeat(vec_str_NAME_OF_CHROMOSOME_OR_SCAFFOLD, inner=n_int_allele_count),
                    #            repeat(vec_int_POSITION, inner=n_int_allele_count),
                    #            mat_int_ALLELE_COUNTS
                    #           )
                    # out = join([join(x,',') for x in eachrow(OUT)], '\n')
                    # file = open(str_filename_output, "a")
                    # write(file, string(out, '\n'))
                    # close(file)
                    fun_writeout_inrun(vec_str_NAME_OF_CHROMOSOME_OR_SCAFFOLD,
                                    vec_int_POSITION,
                                    mat_int_ALLELE_COUNTS,
                                    n_int_allele_count,
                                    str_filename_output)
                end
            end
        end
    end
    close(FILE)
    # return(vec_str_NAME_OF_CHROMOSOME_OR_SCAFFOLD, vec_int_POSITION, mat_int_ALLELE_COUNTS)
    return(0)
end

end