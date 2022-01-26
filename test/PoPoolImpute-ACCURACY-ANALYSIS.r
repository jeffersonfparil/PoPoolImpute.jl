### Analyse imputation accuracies

fname_ctDNA_tar_gz = "PoPoolImpute-ACCURACY-ctDNA.csv.tar.gz" 
system(paste0("tar -xzf ", fname_ctDNA_tar_gz))
fname_ctDNA = paste(rev(tail(rev(unlist(strsplit(fname_ctDNA_tar_gz, "[.]"))),-2)), collapse='.')

dat = read.csv(fname_ctDNA)
str(dat)

fun_filter_data_and_plot = function(data, vec_bool_filter=rep(TRUE, nrow(dat))){
    subdat = data[vec_bool_filter, ]

    RMSE_counts = sqrt(mean((subdat$True_counts - subdat$Imputed_counts)^2))
    RMSE_freqs = sqrt(mean((subdat$True_freqs - subdat$Imputed_freqs)^2))

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

fun_filter_data_and_plot(data=dat)
fun_filter_data_and_plot(data=dat, vec_bool_filter=c((dat$Depth >=100) & (dat$Depth <=1e6)))


### Clean-up
system(paste0("rm ", fname_ctDNA))


