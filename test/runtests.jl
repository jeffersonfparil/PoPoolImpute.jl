using Test
using Pkg

Pkd.add https://github.com/jeffersonfparil/PoPoolImpute.jl.git

using PoPoolImpute

### Uncompress test pileup file
run(`time tar -xvf test.pileup.tar.xz`)

### simulate 10% missing loci in 10% of the pools
run(`time ./simulate_missing_loci_in_pileup_file.sh -f test.pileup -p 0.10 -l 0.10`)


### input and ouput files
str_filename_withMissing = "out_simissing.pileup"
str_filename_output = string("output-imputed-", time(),".syncx")

### impute
@time Test.@test PoPoolImpute.PopPoolImpute(
                                            str_filename_withMissing,
                                            str_filename_output=str_filename_output
                                            )==0

### load imputation output
X = hcat(split.(readlines(str_filename_output), ",")...)
vec_str_NAME_OF_CHROMOSOME_OR_SCAFFOLD = X[1,:]
vec_int_POSITION = parse.(Int, X[2,:])
mat_int_ALLELE_COUNTS = parse.(Int, X[3:end,:])'

### load pileup without missing data, i.e. pileup before simulating missing data
str_filename_noMissing = "test.pileup"
@time _, _, mat_WITHOUT_MISSING = PoPoolImpute.fun_ascii_allele_states_to_counts_per_window(readlines(str_filename_noMissing))

### load input pileup file with the simulated missing data to identify the coordinates of the missing data
@time _, _, mat_WITH_MISSING = PoPoolImpute.fun_ascii_allele_states_to_counts_per_window(readlines(str_filename_withMissing))
mat_bool_missing = ismissing.(mat_WITH_MISSING)

### extract the true and imputed data
@show Y_true = mat_WITHOUT_MISSING[mat_bool_missing]
@show Y_pred = mat_int_ALLELE_COUNTS[mat_bool_missing]
@show length(Y_true)
@show length(Y_pred)

### calculate imputation accuracy
@show n_flt_RMSE = sqrt(sum((Y_true .- Y_pred).^2)/prod(size(Y_true)))

### plot
using UnicodePlots
@show UnicodePlots.scatterplot(Int.(Y_true), Int.(Y_pred), grid=true, color=:white, canvas=BlockCanvas)

### clean-up
rm(str_filename_noMissing)
rm(str_filename_withMissing)
rm(str_filename_output)
