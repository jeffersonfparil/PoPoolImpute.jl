# Analyse imputation accuracies

## Insert model info
```{julia}
using ProgressMeter
# DIR="/data-weedomics-1/ctDNA"
DIR="/data-weedomics-1/Drosophila"
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

## Concatenate and compress
```{sh}
# d=ctDNA
d=Drosophila
DIR=/data-weedomics-1/${d}
cd $DIR
f=$(ls *ACCURACY*.csv | head -n1)
head -n1 $f > PoPoolImpute-OUTPUT-${d}.csv
for f in $(ls *ACCURACY*.csv)
do
    tail -n+2 $f >> PoPoolImpute-OUTPUT-${d}.csv
done
tar -czvf PoPoolImpute-OUTPUT-${d}.csv.tar.gz PoPoolImpute-OUTPUT-${d}.csv
```

## Analyse
```{R}
# d="ctDNA"
d="Drosophila"
fname_tar_gz = paste0("PoPoolImpute-OUTPUT-", d,".csv.tar.gz")
system(paste0("tar -xzf ", fname_tar_gz))
fname = paste(rev(tail(rev(unlist(strsplit(fname_tar_gz, "[.]"))),-2)), collapse='.')

dat = read.csv(fname)
dat$Deviation_counts = abs(dat$True_counts - dat$Imputed_counts)
dat$Deviation_freqs = abs(dat$True_freqs - dat$Imputed_freqs)
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

fun_filter_data_and_plot_deviation = function(data, vec_bool_filter=rep(TRUE, nrow(dat))){
    subdat = data[vec_bool_filter, ]
    par(mfrow=c(2,2))
    plot(subdat$Depth, subdat$Deviation_counts, type="p", pch=19, col=rgb(0.8,0.1,0.2,alpha=0.5))
    grid()
    plot(subdat$Depth, subdat$Deviation_freqs, type="p", pch=19, col=rgb(0.1,0.5,0.1,alpha=0.5))
    grid()
    plot(subdat$True_counts, subdat$Deviation_counts, type="p", pch=19, col=rgb(0.1,0.1,0.6,alpha=0.5))
    grid()
    plot(subdat$True_freqs, subdat$Deviation_freqs, type="p", pch=19, col=rgb(0.9,0.7,0.1,alpha=0.5))
    grid()
}

for (model in unique(dat$Model)){
    for (dist in c("false", "true")){
        if (((model=="Mean") | (model=="OLS")) & (dist=="true")){
            next
        }
        vec_bool_filter=c((dat$Model==model) & (dat$Distance_PCs==dist))
        if (sum(vec_bool_filter)==0){
            next
        }
        png(paste0(model, "-distPCs_", dist, "-scatterplots.png"), width=1000, height=700)
        fun_filter_data_and_plot(data=dat, vec_bool_filter=vec_bool_filter)
        dev.off()
        png(paste0(model, "-distPCs_", dist, "-depth_vs_deviation.png"), width=1000, height=1000)
        fun_filter_data_and_plot_deviation(data=dat, vec_bool_filter=vec_bool_filter)
        dev.off()
    }
}

### Clean-up
system(paste0("rm ", fname))
```

