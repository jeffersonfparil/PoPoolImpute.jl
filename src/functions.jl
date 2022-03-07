using Random
using ProgressMeter
using GLMNet

### FORMATTING CONVENTION
### (1) variable names: snake_case
### (2) structure names: CamelCase
### (3) function names: SCREAMING

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

### Note that this Window struct does not really even need the mutable keyword since its matrix and vector componenets are mutable it seems
mutable struct Window
    chr::Vector{String} ### vector of chromosome names
    pos::Vector{Int}    ### vector of positions
    ref::Vector{Char}   ### vector of reference alleles
    cou::Matrix{Any}    ### n (window size*number of alleles) rows x p (number of pools) columns
    imp::Matrix{Int}    ### number of times a locus has been imputed (corresponds to the elements of cou)
end

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
            if vec_qua_per_pool[k][1] != '*'
                q = Int(vec_qua_per_pool[k][1]) - 33
            else 
                q = 0
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
                    q>=minimum_quality ? D[i] += 1 : nothing
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

# function CLONE(window::Window)::Window
#     Window(copy(window.chr),
#            copy(window.pos),
#            copy(window.ref),
#            copy(window.cou),
#            copy(window.imp))
# end

function IMPUTE!(window::Window; model=["Mean", "OLS", "RR", "LASSO", "GLMNET"][2], distance=true)::Window
    n, p = size(window.cou)
    ### Find the indices of pools with missing data.
    ### These will be used independently and iteratively as our response variables
    idx_pools = sum(ismissing.(window.cou), dims=1)[1,:] .> 0
    ### If we have at least one pool with no missing data, then we proceed with imputation
    if sum(.!idx_pools) >= 1
        ### Impute per pool with missing data
        X = Int.(window.cou[:, .!idx_pools])
        for j in collect(1:p)[idx_pools]
            # j = collect(1:p)[idx_pools][1]
            y = window.cou[:, j]
            idx_loci = ismissing.(y)
            y_train = Int.(y[.!idx_loci])
            X_train = X[.!idx_loci, :]
            nf, pf = size(X_train)
            # ### Distance covariate (only add if the window is within a single chromosome)
            # if distance & (length(unique(window.chr))==1)
            #     pos = window.pos[reshape(.!idx_loci, (7, length(window.pos)))[1,:]]
            #     m = length(pos)
            #     Z = zeros(Int, m, m)
            #     for i in 1:m
            #         for j in 1:m
            #             Z[i,j] = abs(pos[i] - pos[j])
            #         end
            #     end
            #     X_train = hcat(X_train, repeat(Z, inner=(7,1)))
            # end

            ### Train models
            if model == "Mean"
                B = repeat([1/pf], pf)
            elseif model == "OLS"
                B = try
                    X_train \ y_train
                catch
                    missing
                end
            elseif (model == "RR") | (model == "LASSO") | (model == "GLMNET")
                model=="RR" ? alpha=0 : model=="LASSO" ? alpha=1 : alpha=0.5
                B = try
                    try
                        GLMNet.coef(GLMNet.glmnetcv(X, y, alpha=alpha, tol=1e-7)) # equivalend to mod.path.betas[:, argmin(mod)]
                    catch
                        GLMNet.coef(GLMNet.glmnetcv(X, y, alpha=alpha, tol=1e-3)) # equivalend to mod.path.betas[:, argmin(mod)]
                    end
                catch
                    B = missing
                end
            end

            ### Impute
            if !ismissing(B)
                X_valid = X[idx_loci, :]
                y_imputed = Int.(round.(X_valid * B))
                y_imputed[y_imputed .< 0] .= 0 ### collapse negative counts to zero
                window.cou[idx_loci, j] = y_imputed
                window.imp[idx_loci, j] .+= 1
            end
        end
    else
        ### If don't have a single pool with no missing data, then we return the input window without imputing
        nothing
    end
    return(window)
end

### Test
fname = "/home/jeffersonfparil/Documents/PoPoolImpute.jl/test/test-SIMULATED_MISSING.pileup"
window_size = 20
model = "Mean"
distance = false

FILE = open(fname, "r")
i = 0
window = nothing
while !eof(FILE)
    i += 1
    line = PileupLine(i, readline(FILE));
    locus = PARSE(line)
    if window == nothing
        window = [locus]
        while i < window_size
            i += 1
            line = PileupLine(i, readline(FILE));
            locus = PARSE(line)
            push!(window, locus)
        end
        window = PARSE(window)
    end
    IMPUTE!(window, model=model, distance=distance)
    SLIDE!(window, locus=locus)
end
close(FILE)
