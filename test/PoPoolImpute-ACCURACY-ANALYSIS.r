### Analyse imputation accuracies

### Concatenate output and insert model info
```{julia}
using ProgressMeter
DIR="/data-weedomics-1/ctDNA"
cd(DIR)
vec_str_filenames = readdir()
vec_str_filenames = vec_str_filenames[match.(Regex("ACCURACY"), vec_str_filenames) .!= nothing]
@showprogress for str_filename in vec_str_filenames
    # str_filename = vec_str_filenames[1]
    str_dist, str_model, str_rep = split(str_filename, '-')[(end-2):end]
    str_rep = split(str_rep, '.')[1]
    FILE = open(str_filename, "r")
    file = open(string(str_filename, ".temp"), :"w")
    header = readline(FILE)
    header = string("Rep,Model,Distance_PCs,", header, "\n")
    write(file, header)
    while !eof(FILE)
        line = readline(FILE)
        line = string(join([str_rep,str_model,str_dist, line], ','), "\n")
        write(file, line)
    end
    close(FILE)
    close(file)
    mv(string(str_filename, ".temp"), str_filename, force=true)
end
```
```{sh}
f=$(ls *ACCURACY*.csv | head -n1)
head -n1 $f > PoPoolImpute-OUTPUT-ctDNA.csv
for f in $(ls *ACCURACY*.csv)
do
    tail -n+2 $f >> PoPoolImpute-OUTPUT-ctDNA.csv
done
```


fname_ctDNA_tar_gz = "PoPoolImpute-OUTPUT-ctDNA.csv.tar.gz" 
system(paste0("tar -xzf ", fname_ctDNA_tar_gz))
fname_ctDNA = paste(rev(tail(rev(unlist(strsplit(fname_ctDNA_tar_gz, "[.]"))),-2)), collapse='.')

dat = read.csv(fname_ctDNA)
str(dat)

fun_filter_data_and_plot = function(data, vec_bool_filter=rep(TRUE, nrow(dat))){
    subdat = data[vec_bool_filter, ]

    RMSE_counts = sqrt(mean((subdat$True_counts - subdat$Imputed_counts)^2, na.rm=TRUE))
    RMSE_freqs = sqrt(mean((subdat$True_freqs - subdat$Imputed_freqs)^2, na.rm=TRUE))

    par(mfrow=c(1,2))
    plot(subdat$True_counts, subdat$Imputed_counts, type="p", pch=19, col=rgb(0.26,0.63,0.79,alpha=0.5), xlab="Expected", ylab="Imputed", main="Counts")
    grid()
    lines(c(0,40000), c(0,40000), lty=2, lwd=2, col="red")
    legend("topleft", legend=paste0("RMSE = ", round(RMSE_counts)))

    plot(subdat$True_freqs, subdat$Imputed_freqs, type="p", pch=19, col=rgb(0.50,0.80,0.73,alpha=0.5), xlab="Expected", ylab="Imputed", main="Frequencies")
    grid()
    lines(c(0,1), c(0,1), lty=2, lwd=2, col="red")
    legend("topleft", legend=paste0("RMSE = ", round(RMSE_freqs,4)))
}

# fun_filter_data_and_plot(data=dat)
# fun_filter_data_and_plot(data=dat, vec_bool_filter=c((dat$Depth >=100) & (dat$Depth <=1e6)))
for (model in unique(dat$Model)){
    if ((model=="Mean") | (model=="OLS")){
        vec_bool_filter=c((dat$Model==model) & (dat$Depth >=100) & (dat$Depth <=1e6))
        png(paste0(model, "-distPCs_false-scatterplots.png"), width=1000, height=700)
        fun_filter_data_and_plot(data=dat, vec_bool_filter=vec_bool_filter)
        dev.off()
    } else {
        for (dist in c("false", "true")){
            vec_bool_filter=c((dat$Model==model) & (dat$Distance_PCs==dist) & (dat$Depth >=100) & (dat$Depth <=1e6))
            png(paste0(model, "-distPCs_", dist, "-scatterplots.png"), width=1000, height=700)
            fun_filter_data_and_plot(data=dat, vec_bool_filter=vec_bool_filter)
            dev.off()
        }
        
    }
}




dat$Deviation_counts = abs(dat$True_counts - dat$Imputed_counts)
dat$Deviation_freqs = abs(dat$True_freqs - dat$Imputed_freqs)

par(mfrow=c(2,2))
plot(dat$Depth, dat$Deviation_counts, type="p", pch=19, col=rgb(0.8,0.1,0.2,alpha=0.5))
grid()
plot(dat$Depth, dat$Deviation_freqs, type="p", pch=19, col=rgb(0.1,0.5,0.1,alpha=0.5))
grid()
plot(dat$True_counts, dat$Deviation_counts, type="p", pch=19, col=rgb(0.1,0.1,0.6,alpha=0.5))
grid()
plot(dat$True_freqs, dat$Deviation_freqs, type="p", pch=19, col=rgb(0.9,0.7,0.1,alpha=0.5))
grid()


summary(dat$Depth)


### Clean-up
system(paste0("rm ", fname_ctDNA))


