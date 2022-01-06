using UnicodePlots
include("test_1_CI.jl") ### since if one performs: `julia -e include("test/test_2_accuracy.jl")` the current working directory moves automatically to `test/`

P_missing_pools = 0.10
P_missing_loci = 0.10
n_int_number_of_iterations = 10

vec_flt_fraction_missing_imputed = []
vec_flt_RMSE = []

for t in 1:(n_int_number_of_iterations*10)
    while length(vec_flt_fraction_missing_imputed) < n_int_number_of_iterations
        try
            n_flt_fraction_missing_imputed, n_flt_RMSE = fun_sim_impute_check(P_missing_pools=P_missing_pools,
                                                                            P_missing_loci=P_missing_loci,
                                                                            plot=false)
            append!(vec_flt_fraction_missing_imputed, n_flt_fraction_missing_imputed)
            append!(vec_flt_RMSE, n_flt_RMSE)
        catch
            @show "One of the chunks had too many missing data!"
            continue
        end
    end
end

### Plot
sum(vec_flt_fraction_missing_imputed)/length(vec_flt_fraction_missing_imputed)
sum(vec_flt_RMSE)/length(vec_flt_RMSE)

@show UnicodePlots.histogram(Number.(vec_flt_fraction_missing_imputed))
@show UnicodePlots.histogram(Number.(vec_flt_RMSE))

### Clean-up
rm("test.pileup")
