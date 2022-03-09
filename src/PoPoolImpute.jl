module PoPoolImpute

include("functions.jl")
using .functions: SPLIT, IMPUTE
using Dates
using ProgressMeter

function impute(pileup_with_missing::String; window_size::Int=100, model::String=["Mean", "OLS", "RR", "LASSO", "GLMNET"][2], distance::Bool=true, syncx_imputed::String="", threads::Int=2, lines_per_chunk::Int=10_000)::String
    ### Opening remark
    println("")
    println("####################################################################")
    println("PoPoolImpute: imputation of population- and pool-level genotype data")
    println("(version 0.1.0; release 2022/03/10)")
    println("####################################################################")
    println("Input parameters:")
    @show pileup_with_missing
    @show window_size
    @show model
    @show distance
    if threads > 1
        @show threads
        @show lines_per_chunk
    end
    println("####################################################################")
    println(string("Start time: ", Dates.format(now(), "Y-u-dd(E)THH:MM")))
    println("####################################################################")
    ### Define output file if not specified
    if syncx_imputed == ""
        syncx_imputed = string(join(split(pileup_with_missing, '.')[1:(end-1)], '.'), "-IMPUTED.syncx")
    end
    ### Define the full path to the input and output files since calling functions within @distributed loop will revert back to the root directory from where julia was executed from
    if dirname(pileup_with_missing) == ""
        pileup_with_missing = string(pwd(), "/", pileup_with_missing)
    end
    if dirname(syncx_imputed) == ""
        syncx_imputed = string(pwd(), "/", syncx_imputed)
    end
    ### Split input file for parallel processing if we have more than 1 core or thread available
    filenames = [pileup_with_missing]
    if threads > 1
        filenames = SPLIT(pileup_with_missing, lines_per_chunk, window_size)
    end
    ### Impute
    @time filenames_out = @sync @showprogress @distributed (append!) for f in filenames
        filename_imputed = IMPUTE(f, window_size=window_size, model=model, distance=distance)
        [filename_imputed]
    end
    ### Trim-off overhanging windows and merge
    file_out = open(syncx_imputed, "w")
    for i in 1:length(filenames_out)
        if i < length(filenames_out)
            ### trim trailing window from the first and intervening chunks
            lines = 0
            file_in = open(filenames_out[i], "r")
            while !eof(file_in)
                lines += 1
                readline(file_in);
            end
            close(file_in)
            max_line = lines - (window_size*7)
        else
            ### do not trim the last chunk
            max_line = Inf
        end
        file_in = open(filenames_out[i], "r")
        j = 0
        while (!eof(file_in)) & (j < max_line)
            j += 1
            write(file_out, string(readline(file_in), "\n"))
        end
        close(file_in)
        ### clean up
        rm(filenames[i])        ### pileup chunks
        rm(filenames_out[i])    ### syncx chunks
    end
    close(file_out)
    ### Closing remark
    println("")
    println("####################################################################")
    println(string("End time: ", Dates.format(now(), "Y-u-dd(E)THH:MM")))
    println("Imputation successful. Please find the output file:")
    println(syncx_imputed)
    println("####################################################################")
    return(syncx_imputed)
end

end
