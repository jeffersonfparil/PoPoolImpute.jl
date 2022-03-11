### Naming convention:
### (1) variable names: snake_case
### (2) structure names: CamelCase
### (3) function names: SCREAMING

module functions

using GLMNet
using MultivariateStats
using UnicodePlots
using Random
using ProgressMeter

### OBJECTS
struct PileupLine
    index::Int      ### line number
    line::String    ### a line of the pileup file
end

struct LocusAlleleCounts
    chr::String         ### chromosome
    pos::Int            ### position
    ref::Char           ### reference allele
    dep::Vector{Int}    ### vector of depths per pool
    A::Vector{Int}      ### counts per pool of the A allele
    T::Vector{Int}      ### counts per pool of the T allele
    C::Vector{Int}      ### counts per pool of the C allele
    G::Vector{Int}      ### counts per pool of the G allele
    I::Vector{Int}      ### counts per pool of the I allele (insertion)
    D::Vector{Int}      ### counts per pool of the D allele (deletion)
    N::Vector{Int}      ### counts per pool of the N allele (missing)
end

mutable struct Window
    ### Note that this Window struct does not really even need the mutable keyword since its matrix and vector componenets are mutable it seems
    chr::Vector{String} ### vector of chromosome names
    pos::Vector{Int}    ### vector of positions
    ref::Vector{Char}   ### vector of reference alleles
    cou::Matrix{Any}    ### n (window size*number of alleles) rows x p (number of pools) columns
    imp::Matrix{Int}    ### number of times a locus has been imputed (corresponds to the elements of cou)
end

### DATA PARSING AND EXTRACTION
function PARSE(line::PileupLine, minimum_quality=20)::LocusAlleleCounts
    lin = split(line.line, '\t')
    chr = lin[1]
    pos = parse(Int, lin[2])
    ref = lin[3][1]
    dep = parse.(Int, lin[4:3:end])
    sta = lin[5:3:end]
    qua = lin[6:3:end]
    p = length(dep)
    A=[]; T=[]; C=[]; G=[]; I=[]; D=[]; N=[]
    for i in 1:p
        vec_sta_per_pool = split(replace(replace(sta[i], Regex("\\^.") => ""), "\$" => ""), "") ### priot to splitting - remove alignment start and mapping quality strings, as well as alignment end marker
        vec_qua_per_pool = split(qua[i], "")
        j = k = 0
        push!(A, 0); push!(T, 0); push!(C, 0); push!(G, 0); push!(I, 0); push!(D, 0); push!(N, 0)
        while j < length(vec_sta_per_pool)
            j += 1
            k += 1
            # println(string("i=", i))
            # println(string("j=", j))
            # println(string("k=", k))
            
            ### allele read qualities
            q = try
                    if vec_qua_per_pool[k][1] != '*'
                        Int(vec_qua_per_pool[k][1]) - 33
                    else 
                        q = 0
                    end
                catch
                    nothing
                end

            ### allele states
            str_state = uppercase(vec_sta_per_pool[j])
            if (str_state == "+") | (str_state == "-")
                ### remove insertion and deletion sequences
                if str_state == "+"
                    # push!(s, "I")
                    q>=minimum_quality ? I[i] += 1 : nothing
                else
                    # push!(s, 'D')
                    D[i] += 1 ### no need to test for quality because the allele is deleted
                    k -= 1
                end
                l = 1
                count = parse(Int, vec_sta_per_pool[j+l])
                while count != "error"
                    l += 1
                    count = try
                        parse(Int, vec_sta_per_pool[j+l])
                    catch
                        "error"
                    end
                end
                j = j + l + parse(Int, string(vec_sta_per_pool[(j+1):(j+(l-1))]...))
            elseif ((str_state==".") | 
                    (str_state==",") | 
                    (str_state ∈ ["A", "T", "C", "G"])
                   ) & (q>=minimum_quality)
                if str_state ∈ ["A", "T", "C", "G"]
                    a = uppercase(str_state)[1]
                else
                    a = uppercase(ref)[1]
                end
                if a == 'A'
                    A[i] += 1
                elseif a == 'T'
                    T[i] += 1
                elseif a == 'C'
                    C[i] += 1
                elseif a == 'G'
                    G[i] += 1
                else
                    N[i] += 1
                end
            else
                # push!(s, 'N')
                N[i] += 1
            end
            
        end

    end
    return(LocusAlleleCounts(chr, pos, ref, dep, A, T, C, G, I, D, N))
end

function PARSE(window::Vector{LocusAlleleCounts})::Window
    n = length(window)
    p = length(window[1].dep)
    chr = []
    pos = []
    ref = []
    cou = Array{Any,2}(missing, (n*7, p))
    for i in 1:n
        push!(chr, window[i].chr)
        push!(pos, window[i].pos)
        push!(ref, window[i].ref)
        idx = window[i].dep .> 0
        cou[((i-1)*7)+1, idx] = window[i].A[idx]
        cou[((i-1)*7)+2, idx] = window[i].T[idx]
        cou[((i-1)*7)+3, idx] = window[i].C[idx]
        cou[((i-1)*7)+4, idx] = window[i].G[idx]
        cou[((i-1)*7)+5, idx] = window[i].I[idx]
        cou[((i-1)*7)+6, idx] = window[i].D[idx]
        cou[((i-1)*7)+7, idx] = window[i].N[idx]
    end
    imp = zeros((n*7, p))
    return(Window(chr, pos, ref, cou, imp))
end

function EXTRACT(window::Window, locus::Int)::Window
    Window([window.chr[locus]],
           [window.pos[locus]],
           [window.ref[locus]],
           window.cou[(7*(locus-1))+1:(locus*7), :],
           window.imp[(7*(locus-1))+1:(locus*7), :])
end

function EXTRACT(window::Window, loci::UnitRange{Int})::Window
    Window(window.chr[loci],
           window.pos[loci],
           window.ref[loci],
           window.cou[(7*(loci.start-1))+1:(loci.stop*7), :],
           window.imp[(7*(loci.start-1))+1:(loci.stop*7), :])
end

function SLIDE!(window::Window; locus::LocusAlleleCounts)::Window
    new = PARSE([locus])
    window.chr[1:(end-1)] = window.chr[2:end]
    window.pos[1:(end-1)] = window.pos[2:end]
    window.ref[1:(end-1)] = window.ref[2:end]
    window.cou[1:(end-7), :] = window.cou[8:end, :]
    window.imp[1:(end-7), :] = window.imp[8:end, :]
    window.chr[end] = new.chr[1]
    window.pos[end] = new.pos[1]
    window.ref[end] = new.ref[1]
    window.cou[(end-6):end, :] = new.cou[1:7, :]
    window.imp[(end-6):end, :] = zeros(size(new.cou))
    return(window)
end

### IMPUTATE PER WINDOW
function IMPUTE!(window::Window; model::String=["Mean", "OLS", "RR", "LASSO", "GLMNET"][2], distance::Bool=true)::Window
    n, p = size(window.cou)
    ### Find the indices of pools with missing data.
    ### These will be used independently and iteratively as our response variables
    idx_pools = sum(ismissing.(window.cou), dims=1)[1,:] .> 0
    ### If we have at least one pool with no missing data, then we proceed with imputation
    if sum(.!idx_pools) >= 1
        ### Explanatory variables
        X = Int.(window.cou[:, .!idx_pools])

         ### Distance covariate (only add if the window is within a single chromosome)
         if distance & (length(unique(window.chr))==1)
            m = length(window.pos)
            D = zeros(Int, m, m)
            for i in 1:m
                for j in 1:m
                    D[i,j] = abs(window.pos[i] - window.pos[j])
                end
            end
            Z = MultivariateStats.projection(MultivariateStats.fit(PCA, repeat(D, inner=(7,1)); maxoutdim=3)) ### using the first 3 PCs by default
            X = hcat(X, Z)
        end

        for j in collect(1:p)[idx_pools]
            # j = collect(1:p)[idx_pools][1]
            y = window.cou[:, j]
            idx_loci = ismissing.(y)
            y_train = Int.(y[.!idx_loci])
            X_train = X[.!idx_loci, :]
            nf, pf = size(X_train)
           
            ### Train models
            if model == "Mean"
                β = append!([0.0], repeat([1/pf], pf))
            elseif model == "OLS"
                β = try
                    hcat(ones(nf), X_train) \ y_train
                catch
                    try
                        LinearAlgebra.pinv(hcat(ones(nf), X_train)'*hcat(ones(nf), X_train)) * (hcat(ones(nf), X_train)'*y_train)
                    catch
                        missing
                    end
                end
            elseif (model == "RR") | (model == "LASSO") | (model == "GLMNET")
                model=="RR" ? alpha=0 : model=="LASSO" ? alpha=1 : alpha=0.5
                β = try
                    try
                        GLMNet.coef(GLMNet.glmnetcv(hcat(ones(nf), X_train), y_train, alpha=alpha, tol=1e-7)) # equivalend to mod.path.betas[:, argmin(mod)]
                    catch
                        GLMNet.coef(GLMNet.glmnetcv(hcat(ones(nf), X_train), y_train, alpha=alpha, tol=1e-3)) # equivalend to mod.path.betas[:, argmin(mod)]
                    end
                catch
                    β = missing
                end
            end

            ### Impute
            if !ismissing(β)
                X_valid = X[idx_loci, :]
                y_imputed = Int.(round.(hcat(ones(sum(idx_loci)), X_valid) * β))
                # y_imputed[y_imputed .< 0] .= 0 ### collapse negative counts to zero
                y_imputed .-= minimum(y_imputed)

                y_imputed_mean = append!([], ((window.cou[idx_loci, j] .* window.imp[idx_loci, j]) .+ y_imputed) ./ (window.imp[idx_loci, j] .+ 1))
                y_imputed_mean[ismissing.(y_imputed_mean)] = y_imputed

                window.cou[idx_loci, j] = Int.(round.(y_imputed_mean))
                window.imp[idx_loci, j] .+= 1
            end
        end
    else
        ### If don't have a single pool with no missing data, then we return the input window without imputing
        nothing
    end
    return(window)
end

### I/O
function SAVE(window::Window, filename::String)
    OUT = hcat(repeat(window.chr, inner=7),
               repeat(window.pos, inner=7),
               window.cou)
    out = join([join(x,',') for x in eachrow(OUT)], '\n')
    file = open(filename, "a")
    write(file, string(out, '\n'))
    close(file)
end

function PILEUP2SYNCX(filename::String)::String
    file = open(filename, "r")
    filename_output = string(join(split(filename, '.')[1:(end-1)], '.'), ".syncx")
    while !eof(file)
        SAVE(PARSE([PARSE(PileupLine(1, readline(file)))]), filename_output)
    end
    return(filename_output)
end

function SPLIT(filename::String, lines_per_chunk::Int, window_size::Int)::Vector{String}
    file = open(filename, "r")
    c1 = 0; c2 = 0
    i1 = 0; i2 = (lines_per_chunk+window_size+1)

    while !eof(file)
        line = string(readline(file), "\n");
        ### initialise chunk file
        if i1 == 0
            c1 = c2 + 1
            global out1 = open(string(join(split(filename, ".")[1:(end-1)], "."), "-CHUNK_", c1, ".pileup"), "w")
        end
        if i2 == 0
            c2 = c1 + 1
            global out2 = open(string(join(split(filename, ".")[1:(end-1)], "."), "-CHUNK_", c2, ".pileup"), "w")
        end
        ### write into chunk file
        if i1 < (lines_per_chunk+window_size)
            i1 += 1
            write(out1, line)
        end
        if i2 < (lines_per_chunk+window_size)
            i2 += 1
            write(out2, line)
        end
        ### close chunk file
        if i1 == (lines_per_chunk+window_size)
            i1 += 1
            close(out1)
        end
        if i2 == (lines_per_chunk+window_size)
            i2 += 1
            close(out2)
        end
        ### reset chunk size
        if i1 == lines_per_chunk
            i2 = 0
        end
        if i2 == lines_per_chunk
            i1 = 0
        end
    end
    close(file)
    ### Remove tailing files with size less than or equal to the window size
    dir = if dirname(filename) == ""
            "."
        else 
            dirname(filename)
        end
    ls = readdir(dir)
    ls = ls[match.(Regex(join(split(basename(filename), ".")[1:(end-1)], ".")), ls) .!= nothing]
    ls = ls[match.(Regex("CHUNK"), ls) .!= nothing]
    ls = ls[match.(Regex("pileup"), ls) .!= nothing]
    for f in ls
        i = 0
        file = open(f, "r")
        while !eof(file)
            _ = readline(file)
            i += 1
        end
        close(file)
        if i <= window_size
            rm(f)
            ls = ls[ls .!= f]
        end
    end
    ### It is important to have the full path of the chunks since @distributed tasks reverts the working directory to default, i.e. the location where julia was called
    ls = string.(dirname(filename), "/", ls)
    return(ls)
end

### SPARSITY SIMULATION AND CROSS-VALIDATION
function SIMULATESPARSITY(filename; read_length::Int=100, missing_loci_fraction::Float64=0.50, missing_pools_fraction::Float64=0.25, pileup_simulated_missing::String="")
    ###########################################################
    ### TEST
    # filename = "/data-weedomics-1/test_human.pileup"
    # read_length = 150
    # missing_loci_fraction = 0.50
    # missing_pools_fraction = 0.25
    ###########################################################
    if pileup_simulated_missing==""
        pileup_simulated_missing = string(join(split(filename, '.')[1:(end-1)], '.'), "-SIMULATED_MISSING.pileup")
    end
    println("Counting lines.")
    @show loci_count = countlines(filename)
    println("Counting pools.")
    file_temp = open(filename, "r")
    i = -1
    if i <0
        i += 1
        line = readline(file_temp)
    end
    close(file_temp)
    @show pool_count = Int((length(split(line, '\t')) - 3) / 3)
    println("Randomly sampling loci chunks to set to missing.")
    chunk_count = Int(ceil(loci_count / read_length))
    println("Counting the number of loci and pool which will be set to missing.")
    @show missing_loci_count = Int(round(chunk_count*missing_loci_fraction))
    @show missing_pool_count = Int(round(pool_count*missing_pools_fraction))
    ### Proceed if we will be simulating at least 1 missing datapoint
    if missing_loci_count*missing_pool_count > 0
        println("Randomly sample chunks of loci which will be set to missing.")
        random_chunk = sort(randperm(chunk_count)[1:missing_loci_count])
        random_chunk = ((random_chunk .- 1) .* read_length) .+ 1
        println("Open input and output files, and initialise the iterators")
        FILE = open(filename, "r")
        file_out = open(pileup_simulated_missing, "w")
        i = 0
        j = 1
        println("Simulate missing data.")
        pb = ProgressMeter.Progress(loci_count, i)
        while !eof(FILE)
            i += 1; ProgressMeter.next!(pb)
            line = readline(FILE)
            ### while-looping to properly deal with the last read line
            while i == random_chunk[j]
                ### if we reach the last missing chunk, then stop incrementing
                length(random_chunk) == j ? j : j += 1
                ### parse the tab-delimited line
                vec_line = split(line, '\t')
                ### extract the scaffold or chromosome name
                scaffold_or_chromosome = vec_line[1]
                ### randomly choose pools to get missing data
                idx_pool_rand_missing = randperm(pool_count)[1:missing_pool_count]
                idx_pool_rand_missing = (((idx_pool_rand_missing .- 1) .* 3) .+ 1) .+ 3
                position_ini = parse(Int, vec_line[2])
                position_end = position_ini + (read_length - 1) ### less initial position twice since we're counting the initial position as part of the read length and we've already written it before the forst iteration of the while-loop
                while (scaffold_or_chromosome == vec_line[1]) & (parse(Int, vec_line[2]) <= position_end) & !eof(FILE)
                    ### Set to missing each of the randomly sampled pools in the current locus
                    for k in idx_pool_rand_missing
                        vec_line[k:k+2] = ["0", "*", "*"]
                    end
                    ### Write-out the line with simulated missing data
                    line = join(vec_line, '\t')
                    write(file_out, string(line, '\n'))
                    i += 1; ProgressMeter.next!(pb)
                    line = readline(FILE)
                    vec_line = split(line, '\t')
                end
            end
            ### Write-out the line without missing data
            write(file_out, string(line, '\n'))
        end
        close(FILE)
        close(file_out)
        println("##############################################################")
        println("Missing data simulation finished! Please find the output file:")
        println(pileup_simulated_missing)
        println("##############################################################")
        return(pileup_simulated_missing)
    else
        println("Did not simulate any missing data. Because the fraction of missing loci and pools parameter was too low.")
        return(1)
    end
end

function CROSSVALIDATE(syncx_without_missing, syncx_with_missing, syncx_imputed; plot=false, rmse=false, save=false, csv_out="")
    # syncx_without_missing = "test.syncx"
    # syncx_with_missing = "test-SIMULATED_MISSING.syncx"
    # syncx_imputed = "test-IMPUTED.syncx"
    ### NOTE: we should have the same exact locus corresponding per row across these three files
    ### Pools (list of pool indices including the first 2 columns: chr and pos but we need not worry since these won't ever be missing)
    file = open(syncx_without_missing, "r")
    p = collect(1:length(split(readline(file), ',')))
    close(file)
    ### extract expected and imputed allele counts
    file_without_missing = open(syncx_without_missing, "r")
    file_with_missing    = open(syncx_with_missing, "r")
    file_imputed         = open(syncx_imputed, "r")
    expected = []
    imputed = []
    pool = []
    missing_counter = 0
    while !eof(file_without_missing)
        c = split(readline(file_without_missing), ',')
        m = split(readline(file_with_missing), ',')
        i = split(readline(file_imputed), ',')
        missings = (m .== "missing")
        unimputed = (i .== "missing")
        idx = (missings) .& (.!unimputed)
        c = parse.(Int, c[idx])
        i = parse.(Int, i[idx])
        append!(expected, c)
        append!(imputed, i)
        append!(pool, p[idx])
        missing_counter += sum(missings)
    end
    close(file_without_missing)
    close(file_with_missing)
    close(file_imputed)
    ### calculate allele frequencies
    expected_freq = []
    imputed_freq = []
    for i in pool
        idx = pool .== i
        X = reshape(expected[idx], (7, Int(sum(idx)/7)))
        append!(expected_freq, reshape(X ./ sum(X, dims=1), (length(X), )))
        X = reshape(imputed[idx], (7, Int(sum(idx)/7)))
        append!(imputed_freq, reshape(X ./ sum(X, dims=1), (length(X), )))
    end
    ### plot
    if plot
        plot1 = UnicodePlots.scatterplot(Int.(expected), Int.(imputed),
                                         title="Counts",
                                         grid=true, color=:white, canvas=BlockCanvas)
        plot2 = UnicodePlots.scatterplot(Float64.(expected_freq), Float64.(imputed_freq),
                                         title="Frequencies",
                                         grid=true, color=:white, canvas=BlockCanvas)
        @show plot1
        @show plot2
    end
    if length(imputed) == 0
        println(string("No missing data found in ", syncx_with_missing, "."))
        println("Exiting")
        exit()
    end
    ### Root mean square error
    if rmse
        RMSE_count = sqrt(abs(sum((expected .- imputed).^2)/length(expected)))
        RMSE_freq = sqrt(abs(sum((expected_freq .- imputed_freq).^2)/length(expected_freq)))
        imputed_frac = length(imputed) / missing_counter
        println(string("RMSE_count = ", round(RMSE_count, digits=4)))
        println(string("RMSE_freq = ", round(RMSE_freq, digits=4)))
        println(string("Percent imputed = ", round(imputed_frac*100), "%"))
    end
    ### save expected and imputed counts and frequencies
    if save
        if csv_out == ""
            csv_out = string("Imputation_cross_validation_output-", time(), ".csv")
        end
        file_out = open(csv_out, "w")
        for i in 1:length(expected)
            write(file_out, string(join([expected[i], imputed[i], expected_freq[i], imputed_freq[i]], ','), "\n"))
        end
        close(file_out)
    end
    return(expected, imputed, expected_freq, imputed_freq, imputed_frac)
end

function CLONE(window::Window)::Window
    Window(copy(window.chr),
           copy(window.pos),
           copy(window.ref),
           copy(window.cou),
           copy(window.imp))
end

### Main
function IMPUTE(pileup_with_missing::String; window_size::Int=100, model::String=["Mean", "OLS", "RR", "LASSO", "GLMNET"][2], distance::Bool=true, syncx_imputed::String="")::String
    # pileup_with_missing = "/home/jeffersonfparil/Documents/PoPoolImpute.jl/test/test-SIMULATED_MISSING.pileup"
    # window_size = 20
    # model = "LASSO"
    # distance = true
    # syncx_imputed = "test-IMPUTED.syncx"
    ### Output filename
    if syncx_imputed==""
        syncx_imputed = string(join(split(pileup_with_missing, '.')[1:(end-1)], '.'), "-IMPUTED.syncx")
    end
    ### Impute
    file = open(pileup_with_missing, "r")
    i = 0
    window = []
    while !eof(file)
        if window == []
            while i < window_size
                i += 1
                line = PileupLine(i, readline(file));
                locus = PARSE(line)
                push!(window, locus)
            end
            window = PARSE(Array{LocusAlleleCounts}(window))
            IMPUTE!(window, model=model, distance=distance)
            SAVE(EXTRACT(window, 1), syncx_imputed)
        end
        i += 1
        line = PileupLine(i, readline(file));
        locus = PARSE(line)
        SLIDE!(window, locus=locus)
        IMPUTE!(window, model=model, distance=distance)
        SAVE(EXTRACT(window, 1), syncx_imputed)
    end
    close(file)
    SAVE(EXTRACT(window, 2:window_size), syncx_imputed)
    return(syncx_imputed)
end

end

