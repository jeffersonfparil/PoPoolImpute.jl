module functions

using Random
using ProgressMeter
using LinearAlgebra ### Load linear algebra library for the Moore-Penrose pseudoinverse if the automatic solver fails
using MultivariateStats
using Lasso
using GLMNet

### ACSII to allele state parser per locus
function fun_ascii_allele_states_to_counts_per_locus(vec_int_depth, vec_str_allele_state, str_reference_allele, vec_allele_names=["A", "T", "C", "G", "INS", "DEL", "N"])
    int_allele_count = length(vec_allele_names)
    int_pool_count = length(vec_str_allele_state)
    mat_allele_counts_int_pools_x_m_alleles = Array{Any, 2}(missing, int_allele_count, int_pool_count)
    dic_allele_counts = Dict()
    for j in 1:int_pool_count
        # j = 4
        if vec_int_depth[j] > 0
            str_pool_state = replace(vec_str_allele_state[j], Regex("\\^.") => "") ### remove alignment start and mapping quality strings
            str_pool_state = replace(str_pool_state, "\$" => "") ### remove alignment end marker: '$'
            vec_str_split_state = split(str_pool_state, "")
            for allele in vec_allele_names
                dic_allele_counts[allele] = 0
            end
            # println(str_pool_state)
            int_counter = 0
            while int_counter < length(vec_str_split_state)
                int_counter += 1
                # println("########################")
                # println(j)
                # println(int_counter)
                str_state = uppercase(vec_str_split_state[int_counter])
                if (str_state==".") | (str_state==",")
                    dic_allele_counts[uppercase(str_reference_allele)] += 1
                elseif str_state == "+"
                    dic_allele_counts["INS"] += 1
                    n_length_deletion_sequence = parse(Int, vec_str_split_state[int_counter+1])
                    int_counter = int_counter + 1 + 1 + n_length_deletion_sequence ### remove deletion sequence and go to the next state
                elseif str_state == "-"
                    dic_allele_counts["DEL"] += 1
                    n_length_deletion_sequence = parse(Int, vec_str_split_state[int_counter+1])
                    int_counter = int_counter + 1 + n_length_deletion_sequence ### remove deletion sequence
                elseif str_state == "*"
                    dic_allele_counts["DEL"] += 1
                elseif str_state ∈ vec_allele_names
                    dic_allele_counts[str_state] += 1
                else
                    dic_allele_counts["N"] += 1
                end
            end
        else
            for allele in vec_allele_names
                dic_allele_counts[allele] = missing
            end
        end
        for k in 1:length(vec_allele_names)
            # k = 1
            mat_allele_counts_int_pools_x_m_alleles[k, j] = dic_allele_counts[vec_allele_names[k]]
        end
    end
    return(mat_allele_counts_int_pools_x_m_alleles)
end

### Pileup window to allele counts (n pools x (6 loci x m-loci windows))
function fun_ascii_allele_states_to_counts_per_window(vec_str_input, vec_allele_names=["A", "T", "C", "G", "INS", "DEL", "N"])
    int_window_size = length(vec_str_input)
    int_allele_count = length(vec_allele_names)
    int_pool_count = Int((length(split(vec_str_input[1], "\t")) - 3) / 3)
    ### Inititialise matrix of allele counts across pools and loci within the window
    vec_str_name_of_chromosome_or_scaffold = []
    vec_int_position = Int.([])
    mat_int_window_counts = Array{Any,2}(missing, int_allele_count*int_window_size, int_pool_count)
    for i in 1:int_window_size
        # i = 1 ### test 1 specific locus
        # println(i)
        vec_line = split(vec_str_input[i], "\t")
        push!(vec_str_name_of_chromosome_or_scaffold, vec_line[1])
        append!(vec_int_position, parse(Int, vec_line[2]))
        str_reference_allele = vec_line[3]
        ### extract depths, states, and qualities across all pools
        mat_depth_state_quality = reshape(vec_line[4:end], (3, Int(((length(vec_line)-3))/3)))
        vec_int_depth = parse.(Int, mat_depth_state_quality[1,:])
        vec_str_allele_state = mat_depth_state_quality[2,:]
        # vec_str_allele_qualtity = mat_depth_state_quality[3,:]
        n_idx_start = (int_allele_count*(i-1)) + 1
        n_idx_end   = (int_allele_count*(i-0))
        # println("############################")
        mat_int_window_counts[n_idx_start:n_idx_end, :] = fun_ascii_allele_states_to_counts_per_locus(vec_int_depth,
                                                                                                      vec_str_allele_state,
                                                                                                      str_reference_allele,
                                                                                                      vec_allele_names)
    end
    return(vec_str_name_of_chromosome_or_scaffold, vec_int_position, mat_int_window_counts)
end

### Compute pairwise loci distances
function func_pairwise_loci_distances(vec_str_name_of_chromosome_or_scaffold, vec_int_position)
    p = length(vec_int_position)
    Z = zeros(Int, p, p)
    ### Set distance to missing for loci pairs in dfferent chromosomes or scaffolds
    ### This won't result in imputation.
    ### Hence, windows between overlapping chromosomes or scaffolds are not included in the imputation
    ### Maybe this will improve accuracy????!!!
    @inbounds for i in 1:p
        for j in 1:p
            if vec_str_name_of_chromosome_or_scaffold[i] == vec_str_name_of_chromosome_or_scaffold[j]
                Z[i,j] = abs(vec_int_position[i] - vec_int_position[j])
            end
        end
    end
    return(Z)
end

### Multivariate ridge regression which selects the best tuning parameter lambda, and assumes X's first column is ones, i.e for the intercept
function fun_multivariate_ridge_regression(X, Y; flt_ln_lambda_minimum=-5, flt_ln_lambda_maximum=5, int_lambda_count=10, bool_plot=false)
    ### Assumes X's first column is for the intercept
    vec_flt_lambda = exp.(range(flt_ln_lambda_minimum, flt_ln_lambda_maximum, length=int_lambda_count))
    vec_flt_MSE = Float64.([])
    for lambda in vec_flt_lambda
        # lambda = 0.006
        B = MultivariateStats.ridge(Float64.(X), Float64.(Y), lambda, bias=false)
        append!(vec_flt_MSE, sum((X*B - Y).^2)/size(Y,1))
    end
    if bool_plot
        @show UnicodePlots.scatterplot(log.(vec_flt_lambda), vec_flt_MSE)
    end
    MultivariateStats.ridge(Float64.(X),
          Float64.(Y),
          vec_flt_lambda[vec_flt_MSE.==minimum(vec_flt_MSE)][1],
          bias=false)
end

### Impute
function fun_impute_per_window(mat_int_window_counts, vec_str_name_of_chromosome_or_scaffold, vec_int_position, flt_maximum_fraction_of_pools_with_missing=0.5, flt_maximum_fraction_of_loci_with_missing=0.5; bool_use_distance_matrix=false, str_model=["Mean", "OLS", "RR", "LASSO", "GLMNET"][2], flt_glmnet_alpha=0.5)
    ############################################# We're not dealing with depths here because I feel like it is more convenient and provides more flexibility to filter by depth after imputation
    ### TEST
    # using PoPoolImpute
    # cd("/home/jeffersonfparil/Documents/PoPoolImpute.jl/test/")
    # run(`tar -xvf test.pileup.tar.xz`)
    # using Random; Random.seed!(69)
    # PoPoolImpute.functions.fun_simulate_missing("test.pileup",
    #                         n_sequencing_read_length=10,
    #                         flt_maximum_fraction_of_loci_with_missing=0.5,
    #                         flt_maximum_fraction_of_pools_with_missing=0.1,
    #                         str_filename_pileup_simulated_missing="test-SIMULATED_MISSING.pileup")
    # str_filename_input = "test-SIMULATED_MISSING.pileup"
    # int_start_locus = 10
    # int_window_size = 20
    # vec_str_input = readlines(str_filename_input)[int_start_locus:(int_start_locus+int_window_size-1)]
    # vec_allele_names=["A", "T", "C", "G", "INS", "DEL", "N"]
    # @time vec_str_name_of_chromosome_or_scaffold, vec_int_position, mat_int_window_counts = fun_ascii_allele_states_to_counts_per_window(vec_str_input, vec_allele_names)
    # flt_maximum_fraction_of_pools_with_missing = 0.9
    # flt_maximum_fraction_of_loci_with_missing = 0.9
    # bool_use_distance_matrix = false
    # n_bool_window_with_at_least_one_missing_locus = sum(ismissing.(mat_int_window_counts)) > 0
    #############################################

    ### Number of pools
    int_window_size_x_7_alleles, int_pool_count = size(mat_int_window_counts)
    
    ### Find coordinates of missing loci
    mat_bool_missing = ismissing.(mat_int_window_counts)
    
    ### To impute or not to impute
    vec_bool_idx_pools_with_missing_loci = (sum(mat_bool_missing, dims=1) .> 0)[1,:]
    vec_bool_idx_loci_missing = (sum(mat_bool_missing, dims=2) .> 0)[:,1]
    vec_bool_idx_pools_without_missing_loci = .!vec_bool_idx_pools_with_missing_loci
    vec_bool_idx_loci_nomissing = .!vec_bool_idx_loci_missing

    ### Test if we have enough pools with no missing loci to build our imputation model
    n_bool_do_we_have_enough_pools = sum(vec_bool_idx_pools_with_missing_loci) <= (flt_maximum_fraction_of_pools_with_missing*int_pool_count)
    n_bool_do_we_have_enough_loci  = sum(vec_bool_idx_loci_missing)            <= (flt_maximum_fraction_of_loci_with_missing *int_window_size_x_7_alleles)

    if n_bool_do_we_have_enough_pools & n_bool_do_we_have_enough_loci
        ###########################################################################################################
        ### Model the non-missing allele counts in the pools with missing loci
        ### using the allele counts in the same loci in the pools without missing loci;
        ### then use this model (slopes of length or nrow equal to the number of pools without missing loci)
        ### to predict the allele counts of the missing loci in the pools with missing loci. <<<--- phrase this better eh?!?!?!
        ###########################################################################################################
        ### allele counts of pools with missing loci at loci with missing data (m_P non-missing loci x n_M pools with missing loci)
        Y = Number.(mat_int_window_counts[vec_bool_idx_loci_nomissing, vec_bool_idx_pools_with_missing_loci])
        ### allele counts of pools without missing loci at loci with missing data (m_P non-missing loci x n_P pools without missing locus) + first column is for the intercept
        X = Number.(hcat(ones(sum(vec_bool_idx_loci_nomissing)), mat_int_window_counts[vec_bool_idx_loci_nomissing, vec_bool_idx_pools_without_missing_loci]))
        ### predictors of allele counts in the pools with missing loci (n_P pools without missing missing loci x n_M pools with missing loci)
        ### Model the distribution of allele frequencies among the pools with missing data
        ###     as functions of the allele frequencies of the pools without missing data
        bool_overlapping_chromomosomes = length(unique(vec_str_name_of_chromosome_or_scaffold)) > 1
        if bool_use_distance_matrix
            if !bool_overlapping_chromomosomes
                ### Adding distance convariates
                ### Since we are setting distance to missing if the positions are not on the same chromosome of scaffold
                Z = func_pairwise_loci_distances(repeat(vec_str_name_of_chromosome_or_scaffold, inner=7),
                                                 repeat(vec_int_position, inner=7))
                X = Int.(hcat(X[:,1], Z[vec_bool_idx_loci_nomissing,:], X[:,2:end]))
            else
                X = missing ### skip windows with overlapping chromosomes
            end
        end
        if str_model == "Mean"
            B = nothing
        elseif str_model == "OLS"
            B = try
                ### Automatic julia solver (will use qr decomposition for non-square X, i.e. qr(X)\Y)
                X\Y
            catch
                missing
            end
        elseif str_model == "RR"
            B = try
                fun_multivariate_ridge_regression(X, Y)
            catch
                missing
            end
        elseif str_model== "LASSO" ### Note: iterative per column of Y
            B = try
                B = zeros(size(X,2), size(Y,2))
                for j in 1:size(Y,2)
                    B[:,j] = try
                        Lasso.coef(Lasso.fit(LassoModel, Float64.(X[:,2:end]), Float64.(Y[:,j]), standardize=false, intercept=true, cd_tol=1e-7))
                    catch
                        Lasso.coef(Lasso.fit(LassoModel, Float64.(X[:,2:end]), Float64.(Y[:,j]), standardize=false, intercept=true, cd_tol=1e-3))
                    end
                end
                B
            catch
                missing
            end
        elseif str_model== "GLMNET" ### Note: iterative per column of Y
            ### Default: flt_glmnet_alpha=0.5
            try
                B = zeros(size(X,2), size(Y,2))
                for j in 1:size(Y,2)
                    B[:,j] = try
                        GLMNet.coef(GLMNet.glmnetcv(X, Y[:,j], alpha=flt_glmnet_alpha, tol=1e-7)) # equivalend to mod.path.betas[:, argmin(mod)]
                    catch
                        GLMNet.coef(GLMNet.glmnetcv(X, Y[:,j], alpha=flt_glmnet_alpha, tol=1e-3)) # equivalend to mod.path.betas[:, argmin(mod)]
                    end
                end
            catch
                B = missing
            end
        else
            println("Wrong model specified.")
            println("Exiting.")
            exit()
        end
        ### If our solvers fail return missing
        if ismissing(B)
            Y_pred = missing
        elseif str_model == "Mean"
            X_locus_with_missing = mat_int_window_counts[vec_bool_idx_loci_missing, vec_bool_idx_pools_without_missing_loci]
            Y_pred = reshape(repeat(Int.(round.(sum(X_locus_with_missing, dims=2) ./ size(X_locus_with_missing,2))),
                                    outer=sum(vec_bool_idx_pools_with_missing_loci)),
                             (size(X_locus_with_missing,1), sum(vec_bool_idx_pools_with_missing_loci)))

        else
            ### allele counts of pools without missing loci at the loci with with missing data (m_M missing loci x n_P pools without missing loci)
            X_locus_with_missing = hcat(ones(sum(vec_bool_idx_loci_missing)),
                                        mat_int_window_counts[vec_bool_idx_loci_missing, vec_bool_idx_pools_without_missing_loci])
            if bool_use_distance_matrix
                X_locus_with_missing = hcat(X_locus_with_missing[:,1], Z[vec_bool_idx_loci_missing,:], X_locus_with_missing[:,2:end])
            end
            ### prediced allele counts at the loci with missing data (m_M missing loci x n_M pools with missing loci)
            Y_pred = Int.(round.(abs.(X_locus_with_missing * B)))
        end
    else
        Y_pred = missing
    end
    return(Y_pred, vec_bool_idx_pools_with_missing_loci, vec_bool_idx_loci_missing)
    ### Room for improvement by:
    ###     - adding covariates
    ###     - use frequencies and restrict y's to sum up to 1
    ###     - Mixed model
    ###     - Bayesian inference
    ###     - empirical data of Drosophila, Lolium (Arabidopsis?), and human cancer cells pool-seq data
end

### Simple progress bar
function fun_simple_progress_bar(int_current, int_max, int_length=50)
    ######################
    ### TEST
    # int_current = 50
    # int_max = 100
    # int_length = 50
    ######################
    flt_factor = int_length / int_max
    str_progress = repeat("#", Int(round(int_current*flt_factor)))
    str_pending =  repeat("-", Int(round((int_max-int_current)*flt_factor)))
    flt_perc_done = round(int_current * 100 / int_max, digits=2)
    if int_current != int_max
        print("Progress [$str_progress$str_pending] $flt_perc_done% \u001b[1000D")
    else
        ### Keep the progress bar printed out
        print("Progress [$str_progress$str_pending] $flt_perc_done%\n")
    end
end

### split pileup file into chunks
function fun_split_pileup(str_filename_input; int_chunk_count=2, int_chuck_size=1e6, int_window_size=100, n_bool_add_leading_trailing_windows=true)
    ######################
    ### TEST
    # str_filename_input = "test.pileup"
    # int_chunk_count = 2
    # int_chuck_size = 20
    # int_window_size = 10
    # n_bool_add_leading_trailing_windows = true
    ######################
    println(string("Split the input pileup file into: ", int_chunk_count, " chunks of size: ", int_chuck_size, " + 2x", int_window_size, " loci each (at most)."))
    ### Cut up the input file and add leading and/or trailing windows so that we get average imputated frequencies across sliding windows seamlessly
    open(str_filename_input) do FILE
        ### Split into chunks if we are to split the pileup into more than 1 chunk
        if int_chunk_count > 1
            for i in 1:int_chunk_count
                j = 0 ### locus counter per chunk
                ### Open the chunk files (classified as the current, previous, and next chunks in order to append their respective leading and/or trailing windows)
                file_current = open(string(str_filename_input, "-CHUNK_", i), "a")
                if n_bool_add_leading_trailing_windows
                    i > 1                 ? file_previous = open(string(str_filename_input, "-CHUNK_", i-1), "a") : nothing ### first chunk does not have a leading window
                    i < int_chunk_count ? file_next = open(string(str_filename_input, "-CHUNK_", i+1), "a") :     nothing ### last chunk does not have trailing window
                end
                ### Write the line into each chunk file until we reach the chunk size (up to a maximum of chunk size + 2 x window size; accounting for the leading and/or tailing windows)
                while (j < int_chuck_size) & (!eof(FILE))
                    j += 1
                    line = readline(FILE)
                    ### fill up current chunk
                    write(file_current, string(line, '\n'))
                    if n_bool_add_leading_trailing_windows
                        ### add leading window of the current chunk as the trailing window of previous chunk
                        if (j <= int_window_size) & (i > 1)
                            ### Note that the first chunk does not have a leading window
                            write(file_previous, string(line, '\n'))
                        end
                        ### add trailing window of current chunk as the leading window of the next chunk
                        if (j >= (int_chuck_size-(int_window_size-1))) & (i < int_chunk_count)
                            ### Note that the last chunk does not have trailing window
                            write(file_next, string(line, '\n'))
                        end
                    end
                end
                ### close the chunk files
                close(file_current)
                if n_bool_add_leading_trailing_windows
                    i > 1                 ? close(file_previous) : nothing ### first chunk does not have a leading window
                    i < int_chunk_count ? close(file_next) :     nothing ### last chunk does not have trailing window
                end
            end
        else
            println("No splitting required because we have 1 chunk!")
        end
    end
end

### Write file into disk (for writing the output one locus at a time as the imputation algorithm runs iteratively across sliding windows)
function fun_writeout_inrun(vec_str_NAME_OF_CHROMOSOME_OR_SCAFFOLD, vec_int_POSITION, mat_int_ALLELE_COUNTS, int_allele_count, str_filename_output)
    if !isa(vec_str_NAME_OF_CHROMOSOME_OR_SCAFFOLD, Array)
        vec_str_NAME_OF_CHROMOSOME_OR_SCAFFOLD = [vec_str_NAME_OF_CHROMOSOME_OR_SCAFFOLD]
    end
    if !isa(vec_int_POSITION, Array)
        vec_int_POSITION = [vec_int_POSITION]
    end
    OUT = hcat(repeat(vec_str_NAME_OF_CHROMOSOME_OR_SCAFFOLD, inner=int_allele_count),
               repeat(vec_int_POSITION, inner=int_allele_count),
               mat_int_ALLELE_COUNTS)
    out = join([join(x,',') for x in eachrow(OUT)], '\n')
    file = open(str_filename_output, "a")
    write(file, string(out, '\n'))
    close(file)
end

### Single-threaded imputaion
function fun_single_threaded_imputation(str_filename_input; int_window_size=10, flt_maximum_fraction_of_pools_with_missing=0.5, flt_maximum_fraction_of_loci_with_missing=0.5, str_filename_output="output-imputed.syncx", n_bool_skip_leading_window=true, n_bool_skip_trailing_window=true, bool_use_distance_matrix=false, str_model=["Mean", "OLS", "RR", "LASSO", "GLMNET"][2], flt_glmnet_alpha=0.5)
    ###################################################################
    ### TEST
    # cd("/home/jeff/Documents/PoPoolImpute.jl/test")
    # str_filename_input = "out_simissing.pileup"
    # int_window_size = 10
    # flt_maximum_fraction_of_pools_with_missing = 0.5
    # flt_maximum_fraction_of_loci_with_missing = 0.5
    # str_filename_output = "output-imputed.syncx"
    # n_bool_skip_leading_window = true
    # n_bool_skip_trailing_window = true
    ###################################################################
    ### Count the number of loci, alleles, and pools (check if we have the expected number of columns in the first line of the pileup file, i.e. each pool has 3 columns each)
    int_total_loci = countlines(str_filename_input)
    vec_allele_names=["A", "T", "C", "G", "INS", "DEL", "N"]
    int_allele_count = length(vec_allele_names)
    FILE_to_find_pool_count = open(str_filename_input)
    vec_line = split(readline(FILE_to_find_pool_count), "\t")
    close(FILE_to_find_pool_count)
    int_pool_count = (length(vec_line) - 3) / 3
    if int_pool_count == round(int_pool_count)
        int_pool_count = Int(int_pool_count)
    else
        @show "Ooopsss! Pileup file is corrupted."
        @show string("Expected: ", Int(round(int_pool_count)), " pools but got ", int_pool_count, " pools instead.")
        @show "Please check that each pool or population has 3 columns representing the depth, allele state, and allele quality."
    end
    ### Initialise the vectors of scaffold or chromosome names and positions
    vec_str_NAME_OF_CHROMOSOME_OR_SCAFFOLD = []
    vec_int_POSITION = []
    mat_int_ALLELE_COUNTS = nothing
    ### Initialise loci within window counter, vector containing each line (will fit a maximum of int_window_size loci), and open the file
    int_counter_load_first_n_windows_lines_withe_readline = 0
    vec_str_input = []
    FILE = open(str_filename_input)
    ### Iterate per line until we reach the first loci of the last sliding window
    for int_start_locus in 1:((int_total_loci-int_window_size) + 1)
        # int_start_locus = 25
        # @show int_start_locus
        # fun_simple_progress_bar(int_start_locus, (int_total_loci-int_window_size) + 1, 50)
        ### If we already have "int_window_size" lines then just remove the old locus and replace with the next one since we have sliding windows (sliding one locus at a time)
        if int_counter_load_first_n_windows_lines_withe_readline == int_window_size
            vec_str_input[1:(end-1)] = vec_str_input[2:end]
            vec_str_input[end] = readline(FILE) ### read the next line of the file (iterates per line as long as "close(FILE)" has not been executed)
        else
            ### Fill "vec_str_input" so we have "int_window_size" loci
            while int_counter_load_first_n_windows_lines_withe_readline < int_window_size
                push!(vec_str_input, readline(FILE))
                int_counter_load_first_n_windows_lines_withe_readline += 1
            end
        end
        ### Parse each window (contained in "vec_str_input") to extract the chromosome or scaffold names, positions, and the matrix of allele counts with "(int_window_size*vec_allele_names) x int_pool_count" dimensions
        vec_str_name_of_chromosome_or_scaffold, vec_int_position, mat_int_window_counts = fun_ascii_allele_states_to_counts_per_window(vec_str_input, vec_allele_names)
        n_bool_window_with_at_least_one_missing_locus = sum(ismissing.(mat_int_window_counts)) > 0
        ### If we have missing loci then impute, else just add the allele counts without missing information
        if n_bool_window_with_at_least_one_missing_locus
            ### Impute by regressing allele counts of the pools with missing data against the allele counts of the pools without missing data in the window; and then predict the missing allele counts
            mat_imputed, 
            vec_bool_idx_pools_with_missing_loci,
            vec_bool_idx_loci_missing = fun_impute_per_window(mat_int_window_counts,
                                                              vec_str_name_of_chromosome_or_scaffold,
                                                              vec_int_position,
                                                              flt_maximum_fraction_of_pools_with_missing,
                                                              flt_maximum_fraction_of_loci_with_missing,
                                                              bool_use_distance_matrix=bool_use_distance_matrix,
                                                              str_model=str_model,
                                                              flt_glmnet_alpha=flt_glmnet_alpha)
            ### Replace missing data with the imputed allele counts if we were able to impute, i.e. we got at mot most "flt_maximum_fraction_of_pools_with_missing" of the pools with missing loci, and "flt_maximum_fraction_of_loci_with_missing" of the loci with missing data
            if !ismissing(mat_imputed)
                mat_int_window_counts[vec_bool_idx_loci_missing, vec_bool_idx_pools_with_missing_loci] = mat_imputed
            end
        else
            mat_imputed = "No imputation needed since no loci were missing."
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
                int_new_loci_count = sum(vec_bool_idx_loci_new_loci_to_add)
                ### If we have imputed allele counts, then use the average of the imputed allele counts
                if mat_imputed != "No imputation needed since no loci were missing."
                    if int_new_loci_count < int_window_size
                        ### Compute the average imputed allele counts by updating the average given new imputed allele counts
                        vec_idx_bool_loci_missing_less_new_loci = vec_bool_idx_loci_missing[1:(end-(int_allele_count*int_new_loci_count))]
                        mat_int_allele_counts_tail_end_old = mat_int_ALLELE_COUNTS[(end-(int_allele_count*(int_window_size-int_new_loci_count))+1):end, :]
                        mat_int_allele_counts_tail_end_new = mat_int_window_counts[1:(end-(int_allele_count*int_new_loci_count)), :]
                        mat_bool_idx_loci_missing_less_new_locus = reshape(vec_idx_bool_loci_missing_less_new_loci, (int_allele_count, int_window_size-int_new_loci_count))'
                        vec_int_imputed_loci_counter = ones(Int, int_window_size-int_new_loci_count)
                        vec_bool_index_for_vec_int_imputed_loci_counter = (sum(mat_bool_idx_loci_missing_less_new_locus, dims=2) .> 0)[:,1]
                        vec_n = repeat(vec_int_imputed_loci_counter[vec_bool_index_for_vec_int_imputed_loci_counter], inner=int_allele_count)
                        A = mat_int_allele_counts_tail_end_old[vec_idx_bool_loci_missing_less_new_loci, :]
                        B = mat_int_allele_counts_tail_end_new[vec_idx_bool_loci_missing_less_new_loci, :]
                        C = Int.(round.( ( (A.+(B./vec_n)) ./ 2 ) .* ( (2 .* vec_n) ./ (vec_n .+ 1) ) ))
                        ### Update allele counts with the average
                        mat_int_ALLELE_COUNTS[(end-(int_allele_count*(int_window_size-int_new_loci_count))+1):end, :][vec_idx_bool_loci_missing_less_new_loci, :] = C
                        ### Update the counter which we use to compute the updated average allele count
                        vec_int_imputed_loci_counter = vec_int_imputed_loci_counter .+ vec_bool_index_for_vec_int_imputed_loci_counter
                        vec_int_imputed_loci_counter[1:(end-1)] = vec_int_imputed_loci_counter[2:(end-0)]
                        vec_int_imputed_loci_counter[end] = 1
                    end
                end
                ### Save the file per window because it's nice to have the output written into disk rather than memory in case anything unsavory happens prior to finishing the entire job - then at least we'll have a partial output rather than nothing at all
                if (int_start_locus >= 2)
                    ### Save a locus once we're done with the trailing end of the previous window
                    if n_bool_skip_leading_window == false
                        ### save leading window
                        fun_writeout_inrun(vec_str_NAME_OF_CHROMOSOME_OR_SCAFFOLD[1],
                                        vec_int_POSITION[1],
                                        mat_int_ALLELE_COUNTS[1:int_allele_count, :],
                                        int_allele_count,
                                        str_filename_output)
                    else
                        ### do not save the leading window
                        if int_start_locus > (int_window_size+1)
                            fun_writeout_inrun(vec_str_NAME_OF_CHROMOSOME_OR_SCAFFOLD[1],
                                        vec_int_POSITION[1],
                                        mat_int_ALLELE_COUNTS[1:int_allele_count, :],
                                        int_allele_count,
                                        str_filename_output)
                        end
                    end
                end
                ### Update the matrix of allele counts, and vectors of loci coordinates with the new loci keeping the size constant by removing the loci out of the window and adding the new loci entering the window
                if int_new_loci_count < int_window_size
                    ### If the new and old windows overlap as expected with sliding windows, then update the "mat_int_window_counts" according to the extent of the overlap
                    mat_int_ALLELE_COUNTS[1:(end-(int_allele_count*int_new_loci_count)), :] = mat_int_ALLELE_COUNTS[((int_allele_count*int_new_loci_count)+1):end, :]
                    mat_int_ALLELE_COUNTS[((end-(int_allele_count*int_new_loci_count))+1):end, :] = mat_int_window_counts[repeat(vec_bool_idx_loci_new_loci_to_add, inner=int_allele_count), :]
                    vec_str_NAME_OF_CHROMOSOME_OR_SCAFFOLD[1:(end-int_new_loci_count)] = vec_str_NAME_OF_CHROMOSOME_OR_SCAFFOLD[(int_new_loci_count+1):end]
                    vec_str_NAME_OF_CHROMOSOME_OR_SCAFFOLD[((end-int_new_loci_count)+1):end] = vec_str_name_of_chromosome_or_scaffold[vec_bool_idx_loci_new_loci_to_add]
                    vec_int_POSITION[1:(end-int_new_loci_count)] = vec_int_POSITION[(int_new_loci_count+1):end]
                    vec_int_POSITION[((end-int_new_loci_count)+1):end] = vec_int_position[vec_bool_idx_loci_new_loci_to_add]
                else
                    ### If the old and new windows do not overlap, i.e. when "mat_int_window_counts" loci were skipped because of too much missing data which rendered imputation impossible,
                    ### then replace the matrix of allele counts, and vectors of loci coordinates with the new loci
                    mat_int_ALLELE_COUNTS[1:end,:] = mat_int_window_counts
                    vec_str_NAME_OF_CHROMOSOME_OR_SCAFFOLD[1:end] = vec_str_name_of_chromosome_or_scaffold
                    vec_int_POSITION[1:end] = vec_int_position
                end
                ### If we reach the end of the file offset by one window then save the last window
                if int_start_locus == ((int_total_loci-int_window_size) + 1)
                    ### save trailing window
                    if n_bool_skip_trailing_window == false
                        for i in 1:length(vec_str_NAME_OF_CHROMOSOME_OR_SCAFFOLD)
                            fun_writeout_inrun(vec_str_NAME_OF_CHROMOSOME_OR_SCAFFOLD[i],
                                            vec_int_POSITION[i],
                                            mat_int_ALLELE_COUNTS[(((i-1)*int_allele_count)+1):(i*int_allele_count), :],
                                            int_allele_count,
                                            str_filename_output)
                        end
                    end
                end
            end
        end
    end
    ### Note: loci which cannot be imputed are not included in the output file
    ### Close input file
    close(FILE)
    ### Return zero to indicate that all is well
    return(0)
end

###@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
### FUNCTIONS FOR SIMULATING MISSING DATA:

### Filter pileup file: remove loci with at least 1 missing data point
function fun_filter_pileup(str_filename_input; flt_maximum_missing=0.50, str_filename_filtered_pileup=".")
    ######################
    ### TEST
    # str_filename_input = "test.pileup"
    # flt_maximum_missing = 0.50
    ######################
    if str_filename_filtered_pileup == "."
        str_filename_filtered_pileup = string(join(split(str_filename_input, '.')[1:(end-1)], '.'), "-FILTERED_", flt_maximum_missing,".pileup")
    end
    file_filter_pileup = open(str_filename_filtered_pileup, "w")
    int_maximum_missing_threshold = -1
    open(str_filename_input) do FILE
        while !eof(FILE)
            line = readline(FILE)
            vec_counts_quality = split(line, '\t')[4:end]
            if int_maximum_missing_threshold == -1
                int_maximum_missing_threshold = ceil(flt_maximum_missing * (length(vec_counts_quality)/3))
            end
            n_bool_locus_no_depth_in_more_than_N_pools = sum(parse.(Int, vec_counts_quality[collect(1:3:length(vec_counts_quality))]) .== 0) > int_maximum_missing_threshold
            if n_bool_locus_no_depth_in_more_than_N_pools == false
                write(file_filter_pileup, string(line, '\n'))
            end
        end
    end
    close(file_filter_pileup)
end

### simulate missing loci, output a pileup file with missing data, and return the total number of loci
function fun_simulate_missing(str_filename_pileup; n_sequencing_read_length=100, flt_maximum_fraction_of_loci_with_missing=0.50, flt_maximum_fraction_of_pools_with_missing=0.25, str_filename_pileup_simulated_missing=".")
    ###########################################################
    ### TEST
    # str_filename_pileup = "/data-weedomics-1/test_human.pileup"
    # n_sequencing_read_length = 150
    # flt_maximum_fraction_of_loci_with_missing = 0.50
    # flt_maximum_fraction_of_pools_with_missing = 0.25
    ###########################################################
    if str_filename_pileup_simulated_missing=="."
        str_filename_pileup_simulated_missing = string(join(split(str_filename_pileup, '.')[1:(end-1)], '.'), "-SIMULATED_MISSING.pileup")
    end
    println("Counting lines.")
    @show int_loci_count = countlines(str_filename_pileup)
    println("Counting pools.")
    file_temp = open(str_filename_pileup, "r")
    i = -1
    if i <0
        i += 1
        line = readline(file_temp)
    end
    close(file_temp)
    @show int_pool_count = Int((length(split(line, '\t')) - 3) / 3)
    println("Randomly sampling loci chunks to set to missing.")
    int_chunk_count = Int(ceil(int_loci_count / n_sequencing_read_length))
    println("Counting the number of loci and pool which will be set to missing.")
    @show int_missing_loci_count = Int(round(int_chunk_count*flt_maximum_fraction_of_loci_with_missing))
    @show int_missing_pool_count = Int(round(int_pool_count*flt_maximum_fraction_of_pools_with_missing))
    ### Proceed if we will be simulating at least 1 missing datapoint
    if int_missing_loci_count*int_missing_pool_count > 0
        println("Randomly sample chunks of loci which will be set to missing.")
        vec_int_random_chunk = sort(randperm(int_chunk_count)[1:int_missing_loci_count])
        vec_int_random_chunk = ((vec_int_random_chunk .- 1) .* n_sequencing_read_length) .+ 1
        println("Open input and output files, and initialise the interators")
        FILE = open(str_filename_pileup, "r")
        file_out = open(str_filename_pileup_simulated_missing, "w")
        i = 0
        j = 1
        println("Simulate missing data.")
        pb = ProgressMeter.Progress(int_loci_count, i)
        while !eof(FILE)
            i += 1; ProgressMeter.next!(pb)
            line = readline(FILE)
            ### while-looping to properly deal with the last read line
            while i == vec_int_random_chunk[j]
                ### if we reach the last missing chunk, then stop incrementing
                length(vec_int_random_chunk) == j ? j : j += 1
                ### parse the tab-delimited line
                vec_str_line = split(line, '\t')
                ### extract the scaffold or chromosome name
                str_scaffold_or_chromosome = vec_str_line[1]
                ### randomly choose pools to get missing data
                vec_idx_pool_rand_missing = randperm(int_pool_count)[1:int_missing_pool_count]
                vec_idx_pool_rand_missing = (((vec_idx_pool_rand_missing .- 1) .* 3) .+ 1) .+ 3
                int_position_ini = parse(Int, vec_str_line[2])
                int_position_end = int_position_ini + (n_sequencing_read_length - 1) ### less initial position twice since we're counting the initial position as part of the read length and we've already written it before the forst iteration of the while-loop
                while (str_scaffold_or_chromosome == vec_str_line[1]) & (parse(Int, vec_str_line[2]) <= int_position_end) & !eof(FILE)
                    ### Set to missing each of the randomly sampled pools in the current locus
                    for k in vec_idx_pool_rand_missing
                        vec_str_line[k:k+2] = ["0", "*", "*"]
                    end
                    ### Write-out the line with simulated missing data
                    line = join(vec_str_line, '\t')
                    write(file_out, string(line, '\n'))
                    i += 1; ProgressMeter.next!(pb)
                    line = readline(FILE)
                    vec_str_line = split(line, '\t')
                end
            end
            ### Write-out the line without missing data
            write(file_out, string(line, '\n'))
        end
        close(FILE)
        close(file_out)
        println("##############################################################")
        println("Missing data simulation finished! Please find the output file:")
        println(str_filename_pileup_simulated_missing)
        println("##############################################################")
        return(i)
    else
        println("Did not simulate any missing data. Because the fraction of missing loci and pools parameter was too low.")
        return(1)
    end
end

### Find coordinates of missing data in the pileup file
function fun_find_coordinates_of_missing_data(str_filename_withMissing)
    file_with_missing = open(str_filename_withMissing, "r")
    vec_int_idx_missing_loci  = []
    vec_int_idx_missing_pools = []
    vec_str_missing_loci = []
    i = 0
    while !eof(file_with_missing)
        i += 1
        line = readline(file_with_missing)
        vec_str_line = split(line, '\t')
        vec_str_depth_count_quality = vec_str_line[4:end]
        vec_bool_missing = parse.(Int, vec_str_depth_count_quality[1:3:end]) .== 0
        if sum(vec_bool_missing) > 0
            vec_int_idx_pools = findall(vec_bool_missing)
            vec_int_idx_locus = repeat([i], length(vec_int_idx_pools))
            append!(vec_int_idx_missing_loci,  vec_int_idx_locus)
            append!(vec_int_idx_missing_pools, vec_int_idx_pools)
            push!(vec_str_missing_loci, join([join(vec_str_line[1:3], ":"), join(vec_int_idx_pools, ":")], ":"))
        end
    end
    return(vec_int_idx_missing_loci, vec_int_idx_missing_pools, vec_str_missing_loci)
end

### filter output syncx file so it contains only the loci with at least 1 imputed datapoint
function fun_filter_output_syncx(str_filename_output; vec_str_missing_loci, str_filename_syncx_missing_loci_only=".")
    if str_filename_syncx_missing_loci_only == "."
        str_filename_syncx_missing_loci_only = string(join(split(str_filename_output, '.')[1:(end-1)], '.'), "-FILTERED_IMPUTED_LOCI_ONLY.syncx")
    end
    file_imputed = open(str_filename_output, "r")
    file_filtered = open(str_filename_syncx_missing_loci_only, "w")
    i = 1
    vec_str_imputed_loci = []
    println("Counting the number of imputed loci.")
    int_lines_count = countlines(str_filename_output)
    pb = ProgressMeter.Progress(int_lines_count, i)
    println("Filtering imputation output (syncx file) to include only the imputed loci.")
    while !eof(file_imputed)
        ### imputed file loci ID
        line = readline(file_imputed)
        ProgressMeter.next!(pb)
        vec_str_imputed_line = split(line, ',')
        str_scaffold_imputed = vec_str_imputed_line[1]
        int_position_imputed = parse(Int, vec_str_imputed_line[2])
        ### missing loci ID
        vec_str_missing_locus_pool_id = split(vec_str_missing_loci[i], ':')
        str_scaffold_missing = vec_str_missing_locus_pool_id[1]
        int_position_missing = parse(Int, vec_str_missing_locus_pool_id[2])
        ### Assumptions:
        ###     (1) Imputed file and missing loci ID are sorted the same way
        ###     (2) Grouped by scaffold or chromosome
        ###     (3) Position numbers are in increasing order
        ### skip the originally missing loci if it was not imputed (Note: assumes sorted by postion per scaffold or chromosome)
        while (str_scaffold_missing == str_scaffold_imputed) & (int_position_missing < int_position_imputed) & (i < length(vec_str_missing_loci))
            i += 1
            vec_str_missing_locus_pool_id = split(vec_str_missing_loci[i], ':')
            str_scaffold_missing = vec_str_missing_locus_pool_id[1]
            int_position_missing = parse(Int, vec_str_missing_locus_pool_id[2])
        end
        ### Write-out imputed missing loci
        if (str_scaffold_missing == str_scaffold_imputed) & (int_position_missing == int_position_imputed)
            push!(vec_str_imputed_loci, join([str_scaffold_imputed, int_position_imputed], ":"))
            write(file_filtered, string(line, '\n'))
            ### write-out the remaining 6 alleles (A,T,C,G,INS,DEL,N)
            for j in 1:6
                line = readline(file_imputed)
                write(file_filtered, string(line, '\n'))
            end
            i<length(vec_str_missing_loci) ? i += 1 : i = i
        end
    end
    close(file_imputed)
    close(file_filtered)
    return(str_filename_syncx_missing_loci_only, vec_str_imputed_loci)
end

### filter original pileup file to contain only the loci simutated to have missing data and were subsequently imputed
function fun_filter_original_pileup(str_filename_pilelup_no_missing_loci; vec_str_imputed_loci, str_filename_pileup_filtered_imputed_loci=".")
    if str_filename_pileup_filtered_imputed_loci == "."
        str_filename_pileup_filtered_imputed_loci = string(join(split(str_filename_pilelup_no_missing_loci, '.')[1:(end-1)], '.'), "-FILTERED_IMPUTED_LOCI_ONLY.pileup")
    end
    file_orig = open(str_filename_pilelup_no_missing_loci, "r")
    file_imputed = open(str_filename_pileup_filtered_imputed_loci, "w")
    i = 1
    println("Counting the number of loci.")
    int_loci_count = countlines(str_filename_pilelup_no_missing_loci)
    println("Filtering the pileup file to include only the loci which were simulated to be missing and successfully imputed.")
    pb = ProgressMeter.Progress(int_loci_count, i)
    while !eof(file_orig)
        line = readline(file_orig)
        ProgressMeter.next!(pb)
        vec_str_line = split(line, '\t')
        str_scaffold_orig = vec_str_line[1]
        str_position_orig = parse(Int, vec_str_line[2])
        str_scaffold_imputed = split(vec_str_imputed_loci[i], ':')[1]
        str_position_imputed = parse(Int, split(vec_str_imputed_loci[i], ':')[2])
        if str_scaffold_imputed == str_scaffold_orig
            if str_position_imputed == str_position_orig
                i<length(vec_str_imputed_loci) ? i += 1 : i = i
                write(file_imputed, string(line, '\n'))
            end
        end
    end
    close(file_orig)
    close(file_imputed)
    return(str_filename_pileup_filtered_imputed_loci)
end
###@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

end