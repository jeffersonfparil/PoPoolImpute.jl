using Test
using Pkg
using Random
using Distributed
threads = 2 ### github actions virtual machine allocated has only 2 cores
Distributed.addprocs(threads)
# Pkg.add(url="https://github.com/jeffersonfparil/PoPoolImpute.jl.git")
# @everywhere using PoPoolImpute

@everywhere include("/home/jeffersonfparil/Documents/PoPoolImpute.jl/src/PoPoolImpute.jl")

cd("test/")

pileup_without_missing = "test.pileup"
# pileup_with_missing = "test-SIMULATED_MISSING.pileup"
pileup_with_missing = PoPoolImpute.functions.fun_simulate_missing(pileup_without_missing,
                                                                  int_sequencing_read_length=10,
                                                                  flt_maximum_fraction_of_loci_with_missing=0.50,
                                                                  flt_maximum_fraction_of_pools_with_missing=0.25)
for model in ["Mean", "OLS", "RR", "LASSO", "GLMNET"]
    syncx_imputed = PoPoolImpute.impute(pileup_with_missing,
                                        window_size=20,
                                        model=model,
                                        distance=true,
                                        threads=threads,
                                        lines_per_chunk=45)

    @time syncx_without_missing = PoPoolImpute.functions.PILEUP2SYNCX(pileup_without_missing)

    @time syncx_with_missing = PoPoolImpute.functions.PILEUP2SYNCX(pileup_with_missing)

    expected, imputed, expected_freq, imputed_freq = PoPoolImpute.functions.CROSSVALIDATE(syncx_without_missing, syncx_with_missing, syncx_imputed, plot=true, rmse=true, save=true)

    ### clean-up
    for f in [pileup_with_missing, syncx_without_missing, syncx_with_missing, syncx_imputed]
        rm(f)
    end

    ### Remove cross-validation output
    files = readdir()
    for f in files[match.(Regex("Imputation_cross_validation_output-"), files) .!= nothing]
        rm(f)
    end
end