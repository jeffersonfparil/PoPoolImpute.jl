#!/bin/bash
BAM_LIST=$1
REFERENCE_GENOME=$2
MAPQ=$3
BASQ=$4
OUTPUT=$5
####################################################
### TEST:
# DIR=/data-weedomics-1/WeedOmics_database
# find ${DIR}/BAM -name "*.bam" | sort > ${DIR}/MPILEUP/bam_list.txt
# BAM_LIST=${DIR}/MPILEUP/bam_list.txt
# # REFERENCE_GENOME=${DIR}/REFERENCE/lope_V1.0.fasta
# REFERENCE_GENOME=${DIR}/REFERENCE_GENOME/Reference.fasta ### lope_V1.0.fasta fixed to a have single line per scaffold
# MAPQ=0
# BASQ=0
# OUTPUT=${DIR}/MPILEUP/LopeByrneV1-ACC01_62.mpileup
# ###################################################
# cd ${DIR}/MPILEUP
# ./bam_to_mpileup.sh \
#     ${BAM_LIST} \
#     ${REFERENCE_GENOME} \
#     ${MAPQ} \
#     ${BASQ} \
#     ${OUTPUT}
####################################################
time \
samtools mpileup \
    -aa \
    --min-MQ ${MAPQ} \
    --min-BQ ${BASQ} \
    --fasta-ref ${REFERENCE_GENOME} \
    --bam-list ${BAM_LIST} > ${OUTPUT}

### MISC:
# samtools view ACC01.bam | awk '{sum+=$5} END { print "Mean MAPQ =",sum/NR}'
# samtools view ACC01.bam | head -n1 | cut -f11 | awk -l ordchr -v RS='.{1}' '{print ord(RT)+33}'
