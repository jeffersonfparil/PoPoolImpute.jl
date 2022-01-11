#!/bin/bash

DIR_1="/data/weedomics/2.c_60_populations_genotyping"
DIR_2="/data/Lolium/Population_Genetics/2018November_Ryegrass_Collection/VCF"

head -n100 ${DIR_2}/Lolium2019_100X_filtered.mpileup > ${DIR_1}/test.mpileup
grep 'scaffold_31709' ${DIR_2}/Lolium2019.mpileup >> ${DIR_1}/test.mpileup

### Download Pool-seq data from Drosophila and Human cancer cells; and also 

