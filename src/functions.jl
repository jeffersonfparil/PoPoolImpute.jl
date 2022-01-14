module functions

### ACSII to allele state parser per locus
function fun_ascii_allele_states_to_counts_per_locus(vec_int_depth, vec_str_allele_state, str_reference_allele, vec_allele_names=["A", "T", "C", "G", "INS", "DEL", "N"])
    n_int_allele_count = length(vec_allele_names)
    n_int_pool_count = length(vec_str_allele_state)
    mat_allele_counts_n_int_pools_x_m_alleles = Array{Any, 2}(missing, n_int_allele_count, n_int_pool_count)
    dic_allele_counts = Dict()
    for j in 1:n_int_pool_count
        # j = 4
        if vec_int_depth[j] > 0
            str_pool_state = replace(vec_str_allele_state[j], Regex("\\^.") => "") ### remove alignment start and mapping quality strings
            str_pool_state = replace(str_pool_state, "\$" => "") ### remove alignment end marker: '$'
            vec_str_split_state = split(str_pool_state, "")
            for allele in vec_allele_names
                dic_allele_counts[allele] = 0
            end
            # println(str_pool_state)
            n_int_counter = 0
            while n_int_counter < length(vec_str_split_state)
                n_int_counter += 1
                # println("########################")
                # println(j)
                # println(n_int_counter)
                str_state = uppercase(vec_str_split_state[n_int_counter])
                if (str_state==".") | (str_state==",")
                    dic_allele_counts[uppercase(str_reference_allele)] += 1
                elseif str_state == "+"
                    dic_allele_counts["INS"] += 1
                    n_length_deletion_sequence = parse(Int, vec_str_split_state[n_int_counter+1])
                    n_int_counter = n_int_counter + 1 + 1 + n_length_deletion_sequence ### remove deletion sequence and go to the next state
                elseif str_state == "-"
                    dic_allele_counts["DEL"] += 1
                    n_length_deletion_sequence = parse(Int, vec_str_split_state[n_int_counter+1])
                    n_int_counter = n_int_counter + 1 + n_length_deletion_sequence ### remove deletion sequence
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
            mat_allele_counts_n_int_pools_x_m_alleles[k, j] = dic_allele_counts[vec_allele_names[k]]
        end
    end
    return(mat_allele_counts_n_int_pools_x_m_alleles)
end

### Pileup window to allele counts (n pools x (6 loci x m-loci windows))
function fun_ascii_allele_states_to_counts_per_window(vec_str_input, vec_allele_names=["A", "T", "C", "G", "INS", "DEL", "N"])
    n_int_window_size = length(vec_str_input)
    n_int_allele_count = length(vec_allele_names)
    n_int_pool_count = Int((length(split(vec_str_input[1], "\t")) - 3) / 3)
    ### Inititialise matrix of allele counts across pools and loci within the window
    vec_str_name_of_chromosome_or_scaffold = []
    vec_int_position = []
    mat_int_window_counts = Array{Any,2}(missing, n_int_allele_count*n_int_window_size, n_int_pool_count)
    for i in 1:n_int_window_size
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
        n_idx_start = (n_int_allele_count*(i-1)) + 1
        n_idx_end   = (n_int_allele_count*(i-0))
        # println("############################")
        mat_int_window_counts[n_idx_start:n_idx_end, :] = fun_ascii_allele_states_to_counts_per_locus(vec_int_depth,
                                                                                                      vec_str_allele_state,
                                                                                                      str_reference_allele,
                                                                                                      vec_allele_names)
    end
    return(vec_str_name_of_chromosome_or_scaffold, vec_int_position, mat_int_window_counts)
end

### Impute
function fun_impute_per_window(mat_int_window_counts, n_flt_maximum_fraction_of_pools_with_missing=0.5, n_flt_maximum_fraction_of_loci_with_missing=0.5)
    ############################################# We're not dealing with depths here because I feel like it is more convenient to filter by depth after imputation
    ### TEST
    # str_filename_input = "out_simissing.mpileup"
    # n_int_start_locus = 90
    # n_int_window_size = 20
    # vec_str_input = readlines(str_filename_input)[n_int_start_locus:(n_int_start_locus+n_int_window_size-1)]
    # vec_allele_names=["A", "T", "C", "G", "INS", "DEL", "N"]
    # @time vec_str_name_of_chromosome_or_scaffold, vec_int_position, mat_int_window_counts = fun_ascii_allele_states_to_counts_per_window(vec_str_input, vec_allele_names)
    # n_flt_maximum_fraction_of_pools_with_missing = 0.5
    # n_flt_maximum_fraction_of_loci_with_missing = 0.5
    #############################################

    ### Number of pools
    n_int_window_size_x_7_alleles, n_int_pool_count = size(mat_int_window_counts)
    
    ### Find coordinates of missing loci
    mat_bool_missing = ismissing.(mat_int_window_counts)
    
    ### To impute or not to impute
    vec_bool_idx_pools_with_missing_loci = (sum(mat_bool_missing, dims=1) .> 0)[1,:]
    vec_bool_idx_loci_missing = (sum(mat_bool_missing, dims=2) .> 0)[:,1]
    vec_bool_idx_pools_without_missing_loci = .!vec_bool_idx_pools_with_missing_loci
    vec_bool_idx_loci_nomissing = .!vec_bool_idx_loci_missing

    ### Test if we have enough pools with no missing loci to build our imputation model
    n_bool_do_we_have_enough_pools = sum(vec_bool_idx_pools_with_missing_loci) <= (n_flt_maximum_fraction_of_pools_with_missing*n_int_pool_count)
    n_bool_do_we_have_enough_loci  = sum(vec_bool_idx_loci_missing)            <= (n_flt_maximum_fraction_of_loci_with_missing *n_int_window_size_x_7_alleles)

    if n_bool_do_we_have_enough_pools & n_bool_do_we_have_enough_loci
        ###########################################################################################################
        ### Model the non-missing allele counts in the pools with missing loci
        ### using the allele counts in the same loci in the pools without missing loci;
        ### then use this model (slopes of length or nrow equal to the number of pools without missing loci)
        ### to predict the allele counts of the missing loci in the pools with missing loci. <<<--- phrase this better eh?!?!?!
        ###########################################################################################################
        ### allele counts of pools with missing loci at loci with missing data (m_P non-missing loci x n_M pools with missing loci)
        Y = Number.(mat_int_window_counts[vec_bool_idx_loci_nomissing, vec_bool_idx_pools_with_missing_loci])
        ### allele counts of pools without missing loci at loci with missing data (m_P non-missing loci x n_P pools without missing locus)
        X = Number.(mat_int_window_counts[vec_bool_idx_loci_nomissing, vec_bool_idx_pools_without_missing_loci])
        ### predictors of allele counts in the pools with missing loci (n_P pools without missing missing loci x n_M pools with missing loci)
        ### Model the distribution of allele frequencies among the pools with missing data
        ###     as functions of the allele frequencies of the pools without missing data
        B = try
            ### Automatic julia solver
            X\Y
        catch
            ### Moore-Penrose pseudoinverse if the automatic solver fails
            try
                LinearAlgebra.pinv(X'*X)*(X'*Y)
            catch
                missing
            end
        end
        ### If our solver fails return missing
        if ismissing(B)
            Y_pred = missing
        else
            ### allele counts of pools without missing loci at the loci with with missing data (m_M missing loci x n_P pools without missing loci)
            X_locus_with_missing = mat_int_window_counts[vec_bool_idx_loci_missing, vec_bool_idx_pools_without_missing_loci]
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
function fun_simple_progress_bar(n_int_current, n_int_max, n_int_length=50)
    ######################
    ### TEST
    # n_int_current = 50
    # n_int_max = 100
    # n_int_length = 50
    ######################
    flt_factor = n_int_length / n_int_max
    str_progress = repeat("#", Int(round(n_int_current*flt_factor)))
    str_pending =  repeat("-", Int(round((n_int_max-n_int_current)*flt_factor)))
    flt_perc_done = round(n_int_current * 100 / n_int_max, digits=2)
    if n_int_current != n_int_max
        print("Progress [$str_progress$str_pending] $flt_perc_done% \u001b[1000D")
    else
        ### Keep the progress bar printed out
        print("Progress [$str_progress$str_pending] $flt_perc_done%\n")
    end
end

### split pileup file into chunks
function fun_split_pileup(str_filename_input; n_int_chunk_count=2, n_int_chuck_size=1e6, n_int_window_size=100, n_bool_add_leading_trailing_windows=true)
    ######################
    ### TEST
    # str_filename_input = "test.pileup"
    # n_int_chunk_count = 2
    # n_int_chuck_size = 20
    # n_int_window_size = 10
    # n_bool_add_leading_trailing_windows = true
    ######################
    println(string("Split the input pileup file into: ", n_int_chunk_count, " chunks of size: ", n_int_chuck_size, " + 2x", n_int_window_size, " loci each (at most)."))
    ### Cut up the input file and add leading and/or trailing windows so that we get average imputated frequencies across sliding windows seamlessly
    open(str_filename_input) do FILE
        ### Split into chunks if we are to split the pileup into more than 1 chunk
        if n_int_chunk_count > 1
            for i in 1:n_int_chunk_count
                j = 0 ### locus counter per chunk
                ### Open the chunk files (classified as the current, previous, and next chunks in order to append their respective leading and/or trailing windows)
                file_current = open(string(str_filename_input, "-CHUNK_", i), "a")
                if n_bool_add_leading_trailing_windows
                    i > 1                 ? file_previous = open(string(str_filename_input, "-CHUNK_", i-1), "a") : nothing ### first chunk does not have a leading window
                    i < n_int_chunk_count ? file_next = open(string(str_filename_input, "-CHUNK_", i+1), "a") :     nothing ### last chunk does not have trailing window
                end
                ### Write the line into each chunk file until we reach the chunk size (up to a maximum of chunk size + 2 x window size; accounting for the leading and/or tailing windows)
                while (j < n_int_chuck_size) & (!eof(FILE))
                    j += 1
                    line = readline(FILE)
                    ### fill up current chunk
                    write(file_current, string(line, '\n'))
                    if n_bool_add_leading_trailing_windows
                        ### add leading window of the current chunk as the trailing window of previous chunk
                        if (j <= n_int_window_size) & (i > 1)
                            ### Note that the first chunk does not have a leading window
                            write(file_previous, string(line, '\n'))
                        end
                        ### add trailing window of current chunk as the leading window of the next chunk
                        if (j >= (n_int_chuck_size-(n_int_window_size-1))) & (i < n_int_chunk_count)
                            ### Note that the last chunk does not have trailing window
                            write(file_next, string(line, '\n'))
                        end
                    end
                end
                ### close the chunk files
                close(file_current)
                if n_bool_add_leading_trailing_windows
                    i > 1                 ? close(file_previous) : nothing ### first chunk does not have a leading window
                    i < n_int_chunk_count ? close(file_next) :     nothing ### last chunk does not have trailing window
                end
            end
        else
            println("No splitting required because we have 1 chunk!")
        end
    end
end

### Write file into disk (for writing the output one locus at a time as the imputation algorithm runs iteratively across sliding windows)
function fun_writeout_inrun(vec_str_NAME_OF_CHROMOSOME_OR_SCAFFOLD, vec_int_POSITION, mat_int_ALLELE_COUNTS, n_int_allele_count, str_filename_output)
    if !isa(vec_str_NAME_OF_CHROMOSOME_OR_SCAFFOLD, Array)
        vec_str_NAME_OF_CHROMOSOME_OR_SCAFFOLD = [vec_str_NAME_OF_CHROMOSOME_OR_SCAFFOLD]
    end
    if !isa(vec_int_POSITION, Array)
        vec_int_POSITION = [vec_int_POSITION]
    end
    OUT = hcat(repeat(vec_str_NAME_OF_CHROMOSOME_OR_SCAFFOLD, inner=n_int_allele_count),
               repeat(vec_int_POSITION, inner=n_int_allele_count),
               mat_int_ALLELE_COUNTS)
    out = join([join(x,',') for x in eachrow(OUT)], '\n')
    file = open(str_filename_output, "a")
    write(file, string(out, '\n'))
    close(file)
end

### filter pileup file: remove loci with at least 1 missing data point
function fun_filter_pileup(str_filename_input)
    ######################
    ### TEST
    # str_filename_input = "test.pileup"
    ######################
    file_filter_pileup = open(string(str_filename_input, "-FILTERED.pileup"), "w")
    open(str_filename_input) do FILE
        while !eof(FILE)
            line = readline(FILE)
            vec_counts_quality = split(line, '\t')[4:end]
            n_bool_locus_no_depth_in_at_least_1_pool = sum(parse.(Int, vec_counts_quality[collect(1:3:length(vec_counts_quality))]) .== 0) >= 1
            if n_bool_locus_no_depth_in_at_least_1_pool == false
                write(file_filter_pileup, string(line, '\n'))
            end
        end
    end
    close(file_filter_pileup)
end

### Single-threaded impuataion
function fun_single_threaded_imputation(str_filename_input; n_int_window_size=10, n_flt_maximum_fraction_of_pools_with_missing=0.5, n_flt_maximum_fraction_of_loci_with_missing=0.5, str_filename_output="output-imputed.syncx", n_bool_skip_leading_window=true, n_bool_skip_trailing_window=true)
    ###################################################################
    ### TEST
    # cd("/home/jeff/Documents/PoPoolImpute.jl/test")
    # str_filename_input = "out_simissing.pileup"
    # n_int_window_size = 10
    # n_flt_maximum_fraction_of_pools_with_missing = 0.5
    # n_flt_maximum_fraction_of_loci_with_missing = 0.5
    # str_filename_output = "output-imputed.syncx"
    # n_bool_skip_leading_window = true
    # n_bool_skip_trailing_window = true
    ###################################################################
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
        @show "Ooopsss! Pileup file is corrupted."
        @show string("Expected: ", Int(round(n_int_pool_count)), " pools but got ", n_int_pool_count, " pools instead.")
        @show "Please check that each pool or population has 3 columns representing the depth, allele state, and allele quality."
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
        # fun_simple_progress_bar(n_int_start_locus, (n_int_total_loci-n_int_window_size) + 1, 50)
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
                ### If we have imputed allele counts, then use the average of the imputed allele counts
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
                    if n_bool_skip_leading_window == false
                        ### save leading window
                        fun_writeout_inrun(vec_str_NAME_OF_CHROMOSOME_OR_SCAFFOLD[1],
                                        vec_int_POSITION[1],
                                        mat_int_ALLELE_COUNTS[1:n_int_allele_count, :],
                                        n_int_allele_count,
                                        str_filename_output)
                    else
                        ### do not save the leading window
                        if n_int_start_locus > (n_int_window_size+1)
                            fun_writeout_inrun(vec_str_NAME_OF_CHROMOSOME_OR_SCAFFOLD[1],
                                        vec_int_POSITION[1],
                                        mat_int_ALLELE_COUNTS[1:n_int_allele_count, :],
                                        n_int_allele_count,
                                        str_filename_output)
                        end
                    end
                end
                ### Update the matrix of allele counts, and vectors of loci coordinates with the new loci keeping the size constant by removing the loci out of the window and adding the new loci entering the window
                if n_int_new_loci_count < n_int_window_size
                    ### If the new and old windows overlap as expected with sliding windows, then update the "mat_int_window_counts" according to the extent of the overlap
                    mat_int_ALLELE_COUNTS[1:(end-(n_int_allele_count*n_int_new_loci_count)), :] = mat_int_ALLELE_COUNTS[((n_int_allele_count*n_int_new_loci_count)+1):end, :]
                    mat_int_ALLELE_COUNTS[((end-(n_int_allele_count*n_int_new_loci_count))+1):end, :] = mat_int_window_counts[repeat(vec_bool_idx_loci_new_loci_to_add, inner=n_int_allele_count), :]
                    vec_str_NAME_OF_CHROMOSOME_OR_SCAFFOLD[1:(end-n_int_new_loci_count)] = vec_str_NAME_OF_CHROMOSOME_OR_SCAFFOLD[(n_int_new_loci_count+1):end]
                    vec_str_NAME_OF_CHROMOSOME_OR_SCAFFOLD[((end-n_int_new_loci_count)+1):end] = vec_str_name_of_chromosome_or_scaffold[vec_bool_idx_loci_new_loci_to_add]
                    vec_int_POSITION[1:(end-n_int_new_loci_count)] = vec_int_POSITION[(n_int_new_loci_count+1):end]
                    vec_int_POSITION[((end-n_int_new_loci_count)+1):end] = vec_int_position[vec_bool_idx_loci_new_loci_to_add]
                else
                    ### If the old and new windows do not overlap, i.e. when "mat_int_window_counts" loci were skipped because of too much missing data which rendered imputation impossible,
                    ### then replace the matrix of allele counts, and vectors of loci coordinates with the new loci
                    mat_int_ALLELE_COUNTS[1:end,:] = mat_int_window_counts
                    vec_str_NAME_OF_CHROMOSOME_OR_SCAFFOLD[1:end] = vec_str_name_of_chromosome_or_scaffold
                    vec_int_POSITION[1:end] = vec_int_position
                end
                ### If we reach the end of the file offset by one window then save the last window
                if n_int_start_locus == ((n_int_total_loci-n_int_window_size) + 1)
                    ### save trailing window
                    if n_bool_skip_trailing_window == false
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
    end
    ### Note: loci which cannot be imputed are not included in the output file
    ### Close input file
    close(FILE)
    ### Return zero to indicate that all is well
    return(0)
end

end