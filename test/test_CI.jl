using Test
using Pkg
using UnicodePlots
using Distributed
n_int_thread_count = 2 ### guthub actions virtual machine allocated has only 2 cores
Distributed.addprocs(n_int_thread_count)
Pkg.add(url="https://github.com/jeffersonfparil/PoPoolImpute.jl.git")
@everywhere using PoPoolImpute

### Navigate to testing directory
cd("test/")

################################
### TEST LOCALLY: comment-out lines 7 and 8 first
# @everywhere include("/home/jeffersonfparil/Documents/PoPoolImpute.jl/src/PoPoolImpute.jl")
# cd("/home/jeffersonfparil/Documents/PoPoolImpute.jl/test")
################################

### Simulate missing loci, impute, load imputation output, and check imputation accuracy
function fun_sim_impute_check(;P_missing_pools=0.5, P_missing_loci=0.5, n_sequencing_read_length=10, n_int_number_of_iterations=1)
    ### Uncompress test pileup file
    run(`tar -xvf test.pileup.tar.xz`)
    ### Initialise output vectors
    vec_flt_fraction_missing_imputed = []
    vec_flt_RMSE = []
    ### Using a while-try-catch expression to prevent failure when one of the chunks are too sparse due to the stochasticity of "2_simulate_missing_loci_in_pileup_file.sh"
    t = 0
    while t < n_int_number_of_iterations
        ### Simulate 10% missing loci in 10% of the pools
        run(`./3_simulate_missing_loci_in_pileup_file-IMPROVEMENT.sh -f test.pileup -p $P_missing_pools -l $P_missing_loci -r $n_sequencing_read_length`)
        ### Input and ouput files
        str_filename_withMissing = "out_simissing.pileup"
        str_filename_output = string("output-imputed-", time(),".syncx")
        try
            ### Impute
            Test.@test PoPoolImpute.impute(str_filename_withMissing, 
                                                    n_int_window_size=20,
                                                    n_flt_maximum_fraction_of_pools_with_missing=0.5,
                                                    n_flt_maximum_fraction_of_loci_with_missing=0.5,
                                                    str_filename_output=str_filename_output,
                                                    n_int_thread_count=n_int_thread_count)==0
            t += 1
        catch
            @show "At least one of the chunks has too many missing data!"
            ### Clean-up defective chunk/s and the pileup with simulated missing data
            vec_str_fname_chunks = readdir()[match.(r"CHUNK", readdir()) .!= nothing]
            rm.(vec_str_fname_chunks)
            rm(str_filename_withMissing)
            ### Continue with the following iteration of the while-loop
            continue
        end
        ### Load imputation output
        X = hcat(split.(readlines(str_filename_output), ",")...)
        vec_str_NAME_OF_CHROMOSOME_OR_SCAFFOLD = X[1,:]
        vec_int_POSITION = parse.(Int, X[2,:])
        mat_int_ALLELE_COUNTS = parse.(Int, X[3:end,:])'
        ### Load pileup without missing data, i.e. pileup before simulating missing data
        str_filename_noMissing = "test.pileup"
        @time scaf_WITHOUT_MISSING, pos_WITHOUT_MISSING, mat_WITHOUT_MISSING = PoPoolImpute.fun_ascii_allele_states_to_counts_per_window(readlines(str_filename_noMissing))
        ### Load input pileup file with the simulated missing data to identify the coordinates of the missing data
        @time _, _, mat_WITH_MISSING = PoPoolImpute.fun_ascii_allele_states_to_counts_per_window(readlines(str_filename_withMissing))
        ### How much were we able to impute?
        vec_allele_names = ["A", "T", "C", "G", "x_INS", "y_DEL", "z_N"]
        n_int_allele_count = length(vec_allele_names)
        vec_bool_idx_loci_imputation_output_no_repeats = [sum((scaf_WITHOUT_MISSING[i].==vec_str_NAME_OF_CHROMOSOME_OR_SCAFFOLD) .& (pos_WITHOUT_MISSING[i].==vec_int_POSITION))>0 for i in 1:length(scaf_WITHOUT_MISSING)]
        vec_bool_idx_loci_imputation_output = repeat(vec_bool_idx_loci_imputation_output_no_repeats, inner=n_int_allele_count)
        vec_bool_idx_simulated_missing_loci = sum(ismissing.(mat_WITH_MISSING), dims=2).>0
        @show n_flt_fraction_missing_imputed = sum(vec_bool_idx_loci_imputation_output .& vec_bool_idx_simulated_missing_loci) / sum(vec_bool_idx_simulated_missing_loci)
        ### Subset datasets to include only the loci which are included in the imputation output file, i.e. loci with no missing and loci with missing data which we successfully imputed
        mat_WITH_MISSING = mat_WITH_MISSING[vec_bool_idx_loci_imputation_output, :]
        mat_WITHOUT_MISSING = mat_WITHOUT_MISSING[vec_bool_idx_loci_imputation_output, :]
        scaf_WITHOUT_MISSING = scaf_WITHOUT_MISSING[vec_bool_idx_loci_imputation_output_no_repeats]
        pos_WITHOUT_MISSING = pos_WITHOUT_MISSING[vec_bool_idx_loci_imputation_output_no_repeats]
        ### Sort datasets according to loci and alleles
        vec_int_order = sortperm(string.(repeat(scaf_WITHOUT_MISSING, inner=n_int_allele_count), ["-"], repeat(pos_WITHOUT_MISSING, inner=n_int_allele_count), ["-"], repeat(vec_allele_names, outer=length(scaf_WITHOUT_MISSING))))
        mat_WITH_MISSING = mat_WITH_MISSING[vec_int_order, :]
        mat_WITHOUT_MISSING = mat_WITHOUT_MISSING[vec_int_order, :]
        vec_int_order2 = sortperm(string.(vec_str_NAME_OF_CHROMOSOME_OR_SCAFFOLD, ["-"], vec_int_POSITION, ["-"], repeat(vec_allele_names, outer=length(scaf_WITHOUT_MISSING))))
        mat_int_ALLELE_COUNTS = mat_int_ALLELE_COUNTS[vec_int_order2, :]
        ### Identify coordinates of the simulated missing data
        mat_bool_missing = vec(ismissing.(mat_WITH_MISSING))
        ### Extract the true and imputed data
        Y_true = vec(mat_WITHOUT_MISSING)[mat_bool_missing]
        Y_pred = vec(mat_int_ALLELE_COUNTS)[mat_bool_missing]
        if n_int_number_of_iterations == 1
            @show UnicodePlots.scatterplot(Int.(Y_true), Int.(Y_pred), grid=true, color=:white, canvas=BlockCanvas)
        end
        ### Calculate imputation accuracy
        @show n_flt_RMSE = sqrt(sum((Y_true .- Y_pred).^2)/prod(size(Y_true)))
        ### Append fraction of imputed missing data, and RMSE into the the output vectors
        append!(vec_flt_fraction_missing_imputed, n_flt_fraction_missing_imputed)
        append!(vec_flt_RMSE, n_flt_RMSE)
        ### Clean-up
        rm(str_filename_withMissing)
        rm(str_filename_output)
    end
    if n_int_number_of_iterations > 1
        @show sum(vec_flt_fraction_missing_imputed)/length(vec_flt_fraction_missing_imputed)
        @show sum(vec_flt_RMSE)/length(vec_flt_RMSE)
        @show UnicodePlots.histogram(Number.(vec_flt_fraction_missing_imputed))
        @show UnicodePlots.histogram(Number.(vec_flt_RMSE))
    end
    ### Clean-up
    rm("test.pileup")
    ### Output
    return(0)
end


### MISC: USING OTHER DATASETS
### Filter pileups to include no missing data
# vec_str_pileup = ["Drosophila/Drosophila.mpileup",
#                   "Human/Human.mpileup"]
# for str_pileup in vec_str_pileup
#     PoPoolImpute.fun_filter_pileup(str_pileup, flt_maximum_missing=0.5)
#     PoPoolImpute.fun_filter_pileup(str_pileup, flt_maximum_missing=0.0)
# end


### Simulat missing data
julia

using Random
using ProgressMeter

str_filename_pileup = "/data-weedomics-1/test_human.pileup"
n_sequencing_read_length = 150
n_flt_maximum_fraction_of_loci_with_missing = 0.50
n_flt_maximum_fraction_of_pools_with_missing = 0.25

str_filename_pileup_simulated_missing = string(join(split(str_filename_pileup, '.')[1:(end-1)], '.'), "-SIMULATED_MISSING.pileup")


println("Counting lines.")
n_int_loci_count = countlines(str_filename_pileup)

println("Counting pools.")
file_temp = open(str_filename_pileup, "r")
i = -1
if i <0
    i += 1
    line = readline(file_temp)
end
close(file_temp)
n_int_pool_count = Int((length(split(line, '\t')) - 3) / 3)

println("Randomly sampling loci chunks to set to missing.")
n_int_chunk_count = Int(ceil(n_int_loci_count / n_sequencing_read_length))

n_int_missing_loci_count = Int(round(n_int_chunk_count*n_flt_maximum_fraction_of_loci_with_missing))
n_int_missing_pool_count = Int(round(n_int_pool_count*n_flt_maximum_fraction_of_pools_with_missing))


vec_int_random_chunk = sort(randperm(n_int_chunk_count)[1:n_int_missing_loci_count])
vec_int_random_chunk = ((vec_int_random_chunk .- 1) .* n_sequencing_read_length) .+ 1

FILE = open(str_filename_pileup, "r")
file_out = open(str_filename_pileup_simulated_missing, "w")
i = 0
j = 1

pb = ProgressMeter.Progress(n_int_loci_count, i)

while !eof(FILE)
    i += 1; ProgressMeter.next!(pb)
    line = readline(FILE)
    ### while-looping to properly deal with the last read line
    while i == vec_int_random_chunk[j]
        length(vec_int_random_chunk) == j ? j : j += 1
        vec_str_line = split(line, '\t')
        
         ### randomly choose pools to get missing data
         vec_idx_pool_rand_missing = randperm(n_int_pool_count)[1:n_int_missing_pool_count]
         vec_idx_pool_rand_missing = (((vec_idx_pool_rand_missing .- 1) .* 3) .+ 1) .+ 3

        n_int_position_ini = parse(Int, vec_str_line[2])
        n_int_position_end = n_int_position_ini + (n_sequencing_read_length - 1) ### less initial position twice since we're counting the initial position as part of the read length and we've already written it before the forst iteration of the while-loop
        while parse(Int, vec_str_line[2]) <= n_int_position_end

            for k in vec_idx_pool_rand_missing
                vec_str_line[k:k+2] = ["0", "*", "*"]
            end
            line = join(vec_str_line, '\t')
            write(file_out, string(line, '\n'))
            i += 1; ProgressMeter.next!(pb)
            line = readline(FILE)
            vec_str_line = split(line, '\t')
        end

        if i == vec_int_random_chunk[j]
            ### if the current line is one of the randomly chosed loci to be simulated as missing
            continue
        end
    end
    write(file_out, string(line, '\n'))
end
close(FILE)
close(file_out)
