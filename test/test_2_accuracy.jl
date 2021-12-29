include("PoPoolImpute.jl")

### Navigate to testing directory
cd("../test/")

### Uncompress test pileup file
run(`time tar -xvf test.pileup.tar.xz`)

### Simulate 10% missing loci in 10% of the pools
run(`time ./2_simulate_missing_loci_in_pileup_file.sh -f test.pileup -p 0.50 -l 0.50`)


### Input and ouput files
str_filename_withMissing = "out_simissing.pileup"
str_filename_output = string("output-imputed-", time(),".syncx")

### Impute
@time PoPoolImpute.PopPoolImpute(str_filename_withMissing,
                                 n_int_window_size=20,
                                 str_filename_output=str_filename_output)

