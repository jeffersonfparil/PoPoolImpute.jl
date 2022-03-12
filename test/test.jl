using Pkg
using Random
using Distributed
Pkg.add(url="https://github.com/jeffersonfparil/PoPoolImpute.jl.git")
@everywhere using PoPoolImpute

# @everywhere include("/home/jeffersonfparil/Documents/PoPoolImpute.jl/src/PoPoolImpute.jl")

function GITHUB_CI_TEST(n=10, s=42, threads=2)
    ### NOTE: github actions virtual machine allocated has only 2 cores
    Distributed.addprocs(threads)

    cd("test/")
    run(`tar -xvf test.pileup.tar.xz`)
    pileup_without_missing = "test.pileup"
    syncx_without_missing = PoPoolImpute.functions.PILEUP2SYNCX(pileup_without_missing)

    Random.seed!(s)
    random_seeds = abs.(Random.rand(Int, n))
    for i in random_seeds
        println("####################################################################")
        println(i)
        Random.seed!(i)
        pileup_with_missing = PoPoolImpute.functions.SIMULATESPARSITY(pileup_without_missing,
                                                                    read_length=10,
                                                                    missing_loci_fraction=0.50,
                                                                    missing_pools_fraction=0.25)
        syncx_with_missing = PoPoolImpute.functions.PILEUP2SYNCX(pileup_with_missing)

        for model in ["Mean", "OLS", "RR", "LASSO", "GLMNET"]
            syncx_imputed = PoPoolImpute.impute(pileup_with_missing,
                                                window_size=20,
                                                model=model,
                                                distance=true,
                                                threads=threads,
                                                lines_per_chunk=45)

            expected, imputed, expected_freq, imputed_freq, imputed_frac = PoPoolImpute.functions.CROSSVALIDATE(syncx_without_missing, syncx_with_missing, syncx_imputed, plot=true, rmse=true, save=true)
            rm(syncx_imputed)

        end 

        ### Clean-up
        rm(pileup_with_missing)
        files = readdir()
        for f in files[match.(Regex("syncx"), files) .!= nothing]
            rm(f)
        end
    end

    files = readdir()
    for f in files[match.(Regex("Imputation_cross_validation_output-"), files) .!= nothing]
        rm(f)
    end
    rm(pileup_without_missing)
end

function EMPIRICAL_TEST(pileup_without_missing, n=10, s=42, threads=20)
    Distributed.addprocs(threads)
    syncx_without_missing = PoPoolImpute.functions.PILEUP2SYNCX(pileup_without_missing)
    Random.seed!(s)
    random_seeds = abs.(Random.rand(Int, n))
    for i in random_seeds
        println("####################################################################")
        println(i)
        Random.seed!(i)
        pileup_with_missing = PoPoolImpute.functions.SIMULATESPARSITY(pileup_without_missing,
                                                                    read_length=10,
                                                                    missing_loci_fraction=0.50,
                                                                    missing_pools_fraction=0.25)
        syncx_with_missing = PoPoolImpute.functions.PILEUP2SYNCX(pileup_with_missing)

        for model in ["Mean", "OLS", "RR", "LASSO", "GLMNET"]
            syncx_imputed = PoPoolImpute.impute(pileup_with_missing,
                                                window_size=1_000,
                                                model=model,
                                                distance=true,
                                                threads=threads,
                                                lines_per_chunk=10_000)

            expected, imputed, expected_freq, imputed_freq, imputed_frac = PoPoolImpute.functions.CROSSVALIDATE(syncx_without_missing, syncx_with_missing, syncx_imputed, plot=true, rmse=true, save=true)
            rm(syncx_imputed)
        end 

        ### Clean-up
        rm(pileup_with_missing)
        files = readdir()
        for f in files[match.(Regex("syncx"), files) .!= nothing]
            rm(f)
        end
    end
    println("Concatenate these output files:")
    files = readdir()
    for f in files[match.(Regex("Imputation_cross_validation_output-"), files) .!= nothing]
        println(f)
    end
end
