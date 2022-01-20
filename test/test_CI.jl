using Test
using Pkg
using ProgressMeter
using UnicodePlots
using Distributed
n_int_thread_count = 2 ### guthub actions virtual machine allocated has only 2 cores
Distributed.addprocs(n_int_thread_count)
# Pkg.add(url="https://github.com/jeffersonfparil/PoPoolImpute.jl.git")
# @everywhere using PoPoolImpute

### Navigate to testing directory
cd("test/")

################################
### TEST LOCALLY: comment-out lines 7 and 8 first
@everywhere include("/home/jeffersonfparil/Documents/PoPoolImpute.jl/src/PoPoolImpute.jl")
cd("/home/jeffersonfparil/Documents/PoPoolImpute.jl/test")
using Random
Random.seed!(123)
################################

### Main test function:
###     (1) simulate missing loci,
###     (2) impute
###     (3) load imputation output, and
###     (4) check imputation accuracy
function fun_sim_impute_check(input="test.pileup.tar.xz"; window_size=20, P_missing_pools=0.5, P_missing_loci=0.5, n_sequencing_read_length=10, n_int_number_of_iterations=1)
    # ############################
    # ### TEST
    # input="test.pileup.tar.xz"
    # window_size=20
    # P_missing_pools=0.5
    # P_missing_loci=0.5
    # n_sequencing_read_length=10
    # n_int_number_of_iterations=1
    # ############################
    if input == "test.pileup.tar.xz"
        ### Uncompress test pileup file
        run(`tar -xvf test.pileup.tar.xz`)
        str_filename_pilelup_no_missing_loci = "test.pileup"
    else
        str_filename_pilelup_no_missing_loci = input
    end
    ### Initialise output vectors
    vec_flt_fraction_missing_imputed = []
    vec_flt_RMSE = []
    ### Using a while-try-catch expression to prevent failure when one of the chunks are too sparse due to the stochasticity of "2_simulate_missing_loci_in_pileup_file.sh"
    t = 0 ### iteration counter
    q = 0 ### error counter
    q_max = 10 ### maximum number of error
    while (t < n_int_number_of_iterations) & ( q <= q_max)
        ### Simulate 10% missing loci in 10% of the pools
        str_filename_withMissing = string(join(split(str_filename_pilelup_no_missing_loci, '.')[1:(end-1)], '.'), "-SIMULATED_MISSING.pileup")
        PoPoolImpute.functions.fun_simulate_missing(str_filename_pilelup_no_missing_loci,
                             n_sequencing_read_length=n_sequencing_read_length,
                             n_flt_maximum_fraction_of_loci_with_missing=P_missing_loci,
                             n_flt_maximum_fraction_of_pools_with_missing=P_missing_pools,
                             str_filename_pileup_simulated_missing=str_filename_withMissing)
        ### Impute (catch errors when the input file is too sparse, and re-run simulation of missing data)
        str_filename_output = string(join(split(str_filename_pilelup_no_missing_loci, '.')[1:(end-1)], '.'), "-IMPUTED-", time(), ".pileup")
        try
            Test.@test PoPoolImpute.impute(str_filename_withMissing, 
                                                    n_int_window_size=window_size,
                                                    n_flt_maximum_fraction_of_pools_with_missing=P_missing_pools,
                                                    n_flt_maximum_fraction_of_loci_with_missing=P_missing_loci,
                                                    str_filename_output=str_filename_output,
                                                    n_int_thread_count=n_int_thread_count)==0
            t += 1
            q = 0 ### reset error counter
        catch
            @show "At least one of the chunks has too many missing data!"
            ### Clean-up defective chunk/s and the pileup with simulated missing data
            vec_str_fname_chunks = readdir()[match.(r"CHUNK", readdir()) .!= nothing]
            rm.(vec_str_fname_chunks)
            rm(str_filename_withMissing)
            ### Continue with the following iteration of the while-loop
            q += 1
            continue
        end
        println("Find coordinates of missing data.")
        @time vec_int_idx_missing_loci,
              vec_int_idx_missing_pools,
              vec_str_missing_loci = PoPoolImpute.functions.fun_find_coordinates_of_missing_data(str_filename_withMissing)
        println("Filter the output syncx file to include only the loci with imputed missing data.")
        @time str_filename_syncx_missing_loci_only, 
              vec_str_imputed_loci = PoPoolImpute.functions.fun_filter_output_syncx(str_filename_output, vec_str_missing_loci=vec_str_missing_loci)
        ### Proceed if we imputed anything
        println("Filter the original pileup file to include only the loci which were simulated to have missing data and were subsequently imputed.")
        @time str_filename_pileup_filtered_imputed_loci = PoPoolImpute.functions.fun_filter_original_pileup(str_filename_pilelup_no_missing_loci, vec_str_imputed_loci=vec_str_imputed_loci)
        println("Load true frequencies.")
        @time vec_str_scaf_WITHOUT_MISSING,
            vec_int_pos_WITHOUT_MISSING,
            mat_int_ALLELE_COUNTS_NO_MISSING = PoPoolImpute.fun_ascii_allele_states_to_counts_per_window(readlines(str_filename_pileup_filtered_imputed_loci))
        println("Load imputed frequencies.")
        X = hcat(split.(readlines(str_filename_syncx_missing_loci_only), ",")...)
        vec_str_NAME_OF_CHROMOSOME_OR_SCAFFOLD = X[1,:]
        vec_int_POSITION = parse.(Int, X[2,:])
        # mat_int_ALLELE_COUNTS = parse.(Int, X[3:end,:])'
        mat_int_ALLELE_COUNTS = parse.(Float64, X[3:end,:])'

        println("Calculate allele frequencies")
        mat_flt_ALLELE_FREQS = copy(mat_int_ALLELE_COUNTS)
        for i in collect(1:7:size(mat_flt_ALLELE_FREQS,1))
            mat_flt_ALLELE_FREQS[i:(i+6),:] = mat_int_ALLELE_COUNTS[i:(i+6),:] ./ sum(mat_int_ALLELE_COUNTS[i:(i+6),:], dims=1)
        end
        mat_flt_ALLELE_FREQS_NO_MISSING = copy(mat_int_ALLELE_COUNTS_NO_MISSING)
        for i in collect(1:7:size(mat_flt_ALLELE_FREQS_NO_MISSING,1))
            mat_flt_ALLELE_FREQS_NO_MISSING[i:(i+6),:] = mat_int_ALLELE_COUNTS_NO_MISSING[i:(i+6),:] ./ sum(mat_int_ALLELE_COUNTS_NO_MISSING[i:(i+6),:], dims=1)
        end


        println("Fraction of missing data that were successfully imputed.")
        @show n_flt_fraction_missing_imputed = length(vec_str_scaf_WITHOUT_MISSING) / length(vec_str_missing_loci)
        println("Remove zero rows (zero frequency alleles.")
        vec_bool_idx_no_freq_alleles = (sum(mat_int_ALLELE_COUNTS_NO_MISSING, dims=2) .== 0)[:,1]
        mat_int_ALLELE_COUNTS_NO_MISSING = mat_int_ALLELE_COUNTS_NO_MISSING[.!vec_bool_idx_no_freq_alleles, :]
        mat_flt_ALLELE_FREQS_NO_MISSING = mat_flt_ALLELE_FREQS_NO_MISSING[.!vec_bool_idx_no_freq_alleles, :]
        vec_str_NAME_OF_CHROMOSOME_OR_SCAFFOLD = vec_str_NAME_OF_CHROMOSOME_OR_SCAFFOLD[.!vec_bool_idx_no_freq_alleles]
        vec_int_POSITION = vec_int_POSITION[.!vec_bool_idx_no_freq_alleles]
        mat_int_ALLELE_COUNTS = mat_int_ALLELE_COUNTS[.!vec_bool_idx_no_freq_alleles, :]
        mat_flt_ALLELE_FREQS = mat_flt_ALLELE_FREQS[.!vec_bool_idx_no_freq_alleles, :]



        println("Keep only the pools (or columns) which refer to the imputed points.")
        i = 1
        vec_int_counts_true    = []
        vec_int_counts_imputed = []
        vec_flt_freqs_true    = []
        vec_flt_freqs_imputed = []
        @showprogress for locus in vec_str_missing_loci
            # locus = vec_str_missing_loci[1]
            vec_str_locus = split(locus, ":")
            str_scaffold = vec_str_locus[1]
            int_postion = parse(Int, vec_str_locus[2])
            if (str_scaffold==vec_str_NAME_OF_CHROMOSOME_OR_SCAFFOLD[i]) & (int_postion==vec_int_POSITION[i])
                for allele in 1:7
                    i<length(vec_str_NAME_OF_CHROMOSOME_OR_SCAFFOLD) ? i += 1 : i = i
                    vec_int_idx_pools = parse.(Int, vec_str_locus[4:end])
                    append!(vec_int_counts_true,    mat_int_ALLELE_COUNTS_NO_MISSING[i, vec_int_idx_pools])
                    append!(vec_int_counts_imputed, mat_int_ALLELE_COUNTS[i, vec_int_idx_pools])
                    append!(vec_flt_freqs_true,    mat_flt_ALLELE_FREQS_NO_MISSING[i, vec_int_idx_pools])
                    append!(vec_flt_freqs_imputed, mat_flt_ALLELE_FREQS[i, vec_int_idx_pools])
                end
            end
        end
        if n_int_number_of_iterations == 1
            println("Scatter plot.")
            plot1 = UnicodePlots.scatterplot(Float64.(vec_int_counts_true),
                                           Float64.(vec_int_counts_imputed),
                                           title="Counts",
                                           grid=true, color=:white, canvas=BlockCanvas)
            plot2 = UnicodePlots.scatterplot(Float64.(vec_flt_freqs_true),
                                           Float64.(vec_flt_freqs_imputed),
                                           title="Frquencies",
                                           grid=true, color=:white, canvas=BlockCanvas)
            @show plot1
            @show plot2
        end
        println("Calculate imputation accuracy.")
        @show n_flt_RMSE_counts = sqrt(sum((vec_int_counts_true .- vec_int_counts_imputed).^2)/prod(size(vec_int_counts_true)))
        @show n_flt_Rflt_freqs = sqrt(sum((vec_flt_freqs_true .- vec_flt_freqs_imputed).^2)/prod(size(vec_flt_freqs_true)))
        ### Append fraction of imputed missing data, and RMSE into the the output vectors
        append!(vec_flt_fraction_missing_imputed, n_flt_fraction_missing_imputed)
        append!(vec_flt_RMSE, n_flt_RMSE_counts)
        ### Clean-up
        rm(str_filename_withMissing)
        rm(str_filename_output)
        rm(str_filename_syncx_missing_loci_only)
        rm(str_filename_pileup_filtered_imputed_loci)
    end
    ### View statistics
    if (n_int_number_of_iterations > 1) & (length(vec_flt_fraction_missing_imputed)>1)
        @show sum(vec_flt_fraction_missing_imputed)/length(vec_flt_fraction_missing_imputed)
        @show sum(vec_flt_RMSE)/length(vec_flt_RMSE)
        @show UnicodePlots.histogram(Number.(vec_flt_fraction_missing_imputed))
        @show UnicodePlots.histogram(Number.(vec_flt_RMSE))
    end
    ### Clean-up
    if input == "test.pileup.tar.xz"
        rm("test.pileup")
    end
    ### Output
    return(0)
end

### MISC: USING OTHER DATASETS
# ### Accuracy assessment 
# using Test
# using Pkg
# using ProgressMeter
# using UnicodePlots
# using Distributed
# n_int_thread_count = 30
# Distributed.addprocs(n_int_thread_count)
# Pkg.add(url="https://github.com/jeffersonfparil/PoPoolImpute.jl.git")
# @everywhere using PoPoolImpute
# ### NOTE!!!!! Manually load fun_sim_impute_check() from above.
# @time fun_sim_impute_check("/data-weedomics-1/test_human.pileup",
#                            window_size=200,
#                            P_missing_pools=0.5,
#                            P_missing_loci=0.5,
#                            n_sequencing_read_length=100,
#                            n_int_number_of_iterations=1)
