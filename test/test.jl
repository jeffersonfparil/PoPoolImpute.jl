pileup_without_missing = ARGS[1]
githubci = parse(Bool, ARGS[2])
n = parse(Int, ARGS[3])
s = parse(Int, ARGS[4])
threads = parse(Int, ARGS[5])
window_size = parse(Int, ARGS[6])
lines_per_chunk = parse(Int, ARGS[7])
### TEST
# pileup_without_missing="test.pileup"
# githubci=true
# n=10
# s=42
# threads=2
# window_size=20
# lines_per_chunk=45
#
# time julia test/test.jl test.pileup true 1 123 2 20 45
#
# DIR=/data-weedomics-1
# time \
# julia ${DIR}/PoPoolImpute.jl/test/test.jl \
#       ${DIR}/ctDNA/ctDNA.mpileup-FILTERED_0.0.pileup \
#       false \
#       10 \
#       42069 \
#       20 \
#       1000 \
#       10000


using Pkg
using Random
using UnicodePlots
using Distributed
Distributed.addprocs(threads)
Pkg.add(url="https://github.com/jeffersonfparil/PoPoolImpute.jl.git")
@everywhere using PoPoolImpute

# @everywhere include("/home/jeffersonfparil/Documents/PoPoolImpute.jl/src/PoPoolImpute.jl")

if githubci
    cd("test/")
    run(`tar -xvf test.pileup.tar.xz`)
else
    cd(dirname(pileup_without_missing))
end

syncx_without_missing = PoPoolImpute.functions.PILEUP2SYNCX(pileup_without_missing)

Random.seed!(s)
random_seeds = abs.(Random.rand(Int, n))
for i in 1:length(random_seeds)
    println("####################################################################")
    println(random_seeds[i])
    Random.seed!(random_seeds[i])
    pileup_with_missing = PoPoolImpute.functions.SIMULATESPARSITY(pileup_without_missing,
                                                                read_length=10,
                                                                missing_loci_fraction=0.50,
                                                                missing_pools_fraction=0.25)
    syncx_with_missing = PoPoolImpute.functions.PILEUP2SYNCX(pileup_with_missing)

    for model in ["Mean", "OLS", "RR", "LASSO", "GLMNET"]
        syncx_imputed = PoPoolImpute.impute(pileup_with_missing,
                                            window_size=window_size,
                                            model=model,
                                            distance=true,
                                            threads=threads,
                                            lines_per_chunk=lines_per_chunk,
                                            syncx_imputed=string("Imputation_cross_validation_output-", model, "-REP_", i , ".syncx"))
        println("Cross-validating")
        csv_accuracy = PoPoolImpute.functions.CROSSVALIDATE(syncx_without_missing,
                                                            syncx_with_missing,
                                                            syncx_imputed)
        ### plot
        file = open(csv_accuracy)
        expected = []
        imputed = []
        expected_freq = []
        imputed_freq = []
        fraction_of_missing_imputed = 0
        while !eof(file)
            line = split(readline(file), ',')
            if line[2] != "" 
                push!(expected, parse(Float64, line[1]))
                push!(imputed, parse(Float64, line[2]))
                push!(expected_freq, parse(Float64, line[3]))
                push!(imputed_freq, parse(Float64, line[4]))
            else
                fraction_of_missing_imputed = parse(Float64, line[1])
            end
        end
        close(file)
        plot1 = UnicodePlots.scatterplot(Int.(expected), Int.(imputed),
                                            title="Counts",
                                            grid=true, color=:white, canvas=BlockCanvas)
        plot2 = UnicodePlots.scatterplot(Float64.(expected_freq), Float64.(imputed_freq),
                                            title="Frequencies",
                                            grid=true, color=:white, canvas=BlockCanvas)
        @show plot1
        @show plot2
        @show RMSE_count = sqrt(sum((expected .- imputed).^2)/length(expected))
        @show RMSE_freqs = sqrt(sum((expected_freq .- imputed_freq).^2)/length(expected_freq))
        @show fraction_of_missing_imputed
        ### Clean-up
        rm(syncx_imputed)
    end 
    ### Clean-up
    rm(pileup_with_missing)
    rm(syncx_with_missing)
end

### Clean-up
rm(syncx_without_missing)

if githubci
    files = readdir()
    for f in files[match.(Regex("Imputation_cross_validation_output-"), files) .!= nothing]
        rm(f)
    end
    rm(pileup_without_missing)
end
