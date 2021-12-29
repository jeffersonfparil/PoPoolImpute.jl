using Test
using Pkg
Pkg.add(url="https://github.com/jeffersonfparil/PoPoolImpute.jl.git")
using PoPoolImpute

### Navigate to testing directory
cd("test/")

### Uncompress test pileup file
run(`time tar -xvf test.pileup.tar.xz`)

### Simulate 10% missing loci in 10% of the pools
run(`time ./simulate_missing_loci_in_pileup_file.sh -f test.pileup -p 0.50 -l 0.50`)


### Input and ouput files
str_filename_withMissing = "out_simissing.pileup"
str_filename_output = string("output-imputed-", time(),".syncx")

### Impute
@time Test.@test PoPoolImpute.PopPoolImpute(
                                            str_filename_withMissing,
                                            n_int_window_size=20,
                                            str_filename_output=str_filename_output
                                            )==0

### Load imputation output
X = hcat(split.(readlines(str_filename_output), ",")...)
vec_str_NAME_OF_CHROMOSOME_OR_SCAFFOLD = X[1,:]
vec_int_POSITION = parse.(Int, X[2,:])
mat_int_ALLELE_COUNTS = parse.(Int, X[3:end,:])'

### Load pileup without missing data, i.e. pileup before simulating missing data
str_filename_noMissing = "test.pileup"
@time scaf_WITHOUT_MISSING, pos_WITHOUT_MISSING, mat_WITHOUT_MISSING = PoPoolImpute.fun_ascii_allele_states_to_counts_per_window(readlines(str_filename_noMissing))

### Load input pileup file with the simulated missing data to identify the coordinates of the missing data
@time _, _, mat_WITH_MISSING = PoPoolImpute.fun_ascii_allele_states_to_counts_per_window(readlines(str_filename_withMissing))

### How much were we able to impute?
vec_allele_names = ["A", "T", "C", "G", "x_INS", "y_DEL", "z_N"]
n_int_allele_count = length(vec_allele_names)
vec_bool_idx_loci_imputation_output_no_repeats = [sum((scaf_WITHOUT_MISSING[i].==vec_str_NAME_OF_CHROMOSOME_OR_SCAFFOLD) .& (pos_WITHOUT_MISSING[i].==vec_int_POSITION))>0 for i in 1:length(scaf_WITHOUT_MISSING)]
vec_bool_idx_loci_imputation_output = repeat(vec_bool_idx_loci_imputation_output_no_repeats, inner=n_int_allele_count)
vec_bool_idx_simulated_missing_loci = sum(ismissing.(mat_WITH_MISSING), dims=2).>0
@show n_flt_fraction_missing_imputed = sum(vec_bool_idx_loci_imputation_output .& vec_bool_idx_simulated_missing_loci) / sum(vec_bool_idx_simulated_missing_loci)

### Subset datasets to include only the loci which are included in the imputation output file, i.e. loci with no missing and loci with missing data which we successfully imputed
mat_WITH_MISSING = mat_WITH_MISSING[vec_bool_idx_loci_imputation_output, :]
mat_WITHOUT_MISSING = mat_WITHOUT_MISSING[vec_bool_idx_loci_imputation_output, :]
scaf_WITHOUT_MISSING = scaf_WITHOUT_MISSING[vec_bool_idx_loci_imputation_output_no_repeats]
pos_WITHOUT_MISSING = pos_WITHOUT_MISSING[vec_bool_idx_loci_imputation_output_no_repeats]

### Sort datasets according to loci and alleles
vec_int_order = sortperm(string.(repeat(scaf_WITHOUT_MISSING, inner=n_int_allele_count), ["-"], repeat(pos_WITHOUT_MISSING, inner=n_int_allele_count), ["-"], repeat(vec_allele_names, outer=length(scaf_WITHOUT_MISSING))))
mat_WITH_MISSING = mat_WITH_MISSING[vec_int_order, :]
mat_WITHOUT_MISSING = mat_WITHOUT_MISSING[vec_int_order, :]
vec_int_order2 = sortperm(string.(vec_str_NAME_OF_CHROMOSOME_OR_SCAFFOLD, ["-"], vec_int_POSITION, ["-"], repeat(vec_allele_names, outer=length(scaf_WITHOUT_MISSING))))
mat_int_ALLELE_COUNTS = mat_int_ALLELE_COUNTS[vec_int_order2, :]

### Identify coordinates of the simulated missing data
mat_bool_missing = ismissing.(mat_WITH_MISSING)

### Extract the true and imputed data
@show Y_true = mat_WITHOUT_MISSING[mat_bool_missing]
@show Y_pred = mat_int_ALLELE_COUNTS[mat_bool_missing]
@show length(Y_true)
@show length(Y_pred)

### Calculate imputation accuracy
@show n_flt_RMSE = sqrt(sum((Y_true .- Y_pred).^2)/prod(size(Y_true)))

### Plot
using UnicodePlots
@show UnicodePlots.scatterplot(Int.(Y_true), Int.(Y_pred), grid=true, color=:white, canvas=BlockCanvas)

### Clean-up
rm(str_filename_noMissing)
rm(str_filename_withMissing)
rm(str_filename_output)
