#!/bin/bash

DIR_1="/data/weedomics/2.c_60_populations_genotyping"
DIR_2="/data/Lolium/Population_Genetics/2018November_Ryegrass_Collection/VCF"

head -n100 ${DIR_2}/Lolium2019_100X_filtered.mpileup > ${DIR_1}/test.mpileup
grep 'scaffold_31709' ${DIR_2}/Lolium2019.mpileup >> ${DIR_1}/test.mpileup

### Download Pool-seq data from Drosophila and Human cancer cells; and also 

sudo apt install -y curlftpfs

### Drosophila
for d in Drosophila-ERR173 Drosophila-ERR115 Human-SRR350
do
    # d=Drosophila-ERR173
    mkdir ${d}
    curlftpfs -r ftp.sra.ebi.ac.uk/vol1/fastq/ERR173 Drosophila-ERR173/
    time find Drosophila-ERR173/ -path '*.fastq.gz' > fname_sequences-${d}.temp
    fusermount -u ${d}
    mv fname_sequences-${d}.temp ${d}
    cd ${d}
    sed -i "s/${d%-*}-/ftp.sra.ebi.ac.uk\/vol1\/fastq\//g" fname_sequences-${d}.temp
    # f=Drosophila-ERR173/ERR173999/ERR173999.fastq.gz
    time parallel wget {} ::: $(cat fname_sequences-${d}.temp)
done

