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
                else
                    dic_allele_counts[str_state] += 1
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
        B = X\Y
        ### allele counts of pools without missing loci at the loci with with missing data (m_M missing loci x n_P pools without missing loci)
        X_locus_with_missing = mat_int_window_counts[vec_bool_idx_loci_missing, vec_bool_idx_pools_without_missing_loci]
        ### prediced allele counts at the loci with missing data (m_M missing loci x n_M pools with missing loci)
        Y_pred = Int.(round.(abs.(X_locus_with_missing * B)))
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

### simple progress bar
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

### In-run wite file into disk
function fun_writeout_inrun(vec_str_NAME_OF_CHROMOSOME_OR_SCAFFOLD, vec_int_POSITION, mat_int_ALLELE_COUNTS, n_int_allele_count, str_filename_output)
    if !isa(vec_str_NAME_OF_CHROMOSOME_OR_SCAFFOLD, Array)
        vec_str_NAME_OF_CHROMOSOME_OR_SCAFFOLD = [vec_str_NAME_OF_CHROMOSOME_OR_SCAFFOLD]
    end
    if !isa(vec_int_POSITION, Array)
        vec_int_POSITION = [vec_int_POSITION]
    end
    OUT = hcat(repeat(vec_str_NAME_OF_CHROMOSOME_OR_SCAFFOLD, inner=n_int_allele_count),
               repeat(vec_int_POSITION, inner=n_int_allele_count),
               mat_int_ALLELE_COUNTS
              )
    out = join([join(x,',') for x in eachrow(OUT)], '\n')
    file = open(str_filename_output, "a")
    write(file, string(out, '\n'))
    close(file)
end

end