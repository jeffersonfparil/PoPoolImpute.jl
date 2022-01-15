#!/bin/bash

echo "########################################"
echo "###                                  ###"
echo "### Using local Lolium Pool-seq data ###"
echo "###                                  ###"
echo "########################################"
DIR_1="/data/weedomics/2.c_60_populations_genotyping"
DIR_2="/data/Lolium/Population_Genetics/2018November_Ryegrass_Collection/VCF"
head -n100 ${DIR_2}/Lolium2019_100X_filtered.mpileup > ${DIR_1}/test.mpileup
grep 'scaffold_31709' ${DIR_2}/Lolium2019.mpileup >> ${DIR_1}/test.mpileup



echo "##################################"
echo "###                            ###"
echo "### Download raw Pool-seq data ###"
echo "###                            ###"
echo "##################################"
### Download Pool-seq data from Drosophila and Human cancer cells; and also 
### NCBI's SRA toolkit and GNU's parallel
sudo apt install -y sra-toolkit parallel
### Configure sra-toolkit to set the temporary folder to be in a bigger directory than root, e.g. a virtual storage drive
mkdir ncbi-sra-tools-cache/
vdb-config -i ### invoke the interactive mode
### (1) Go to CACHE
### (2) Choose the location of user-repository as: `ncbi-sra-tools-cache/`
### (3) Save and exit

echo "#####################################"
echo "### Download Drosophila Pool-seq data"
### Note: 100 bp single-end reads
mkdir Drosophila/ ### https://doi.org/10.1534/genetics.116.197335
### Note: fasterq-dump is the new software set to replace fastq-dump, but in my opinion it is still immature as it lacks some of the features that fastq-dump has, e.g. --gzip
time parallel fastq-dump \
                --gzip \
                --skip-technical \
                --readids \
                --read-filter pass \
                --dumpbase \
                --clip \
                --outdir Drosophila/ \
                {} :::  SRR4478513 \
                        SRR4478514 \
                        SRR4478515 \
                        SRR4478516 \
                        SRR4478517 \
                        SRR4478518 \
                        SRR4478519 \
                        SRR4478520

echo "##########################################"
echo "### Download Human caner Pool-seq data ###"
### Note: 100 bp paired-end exome sequencing reads
mkdir Human/ ### https://www.mdpi.com/2075-1729/12/1/41/htm
time parallel fastq-dump \
                --gzip \
                --skip-technical \
                --readids \
                --read-filter pass \
                --dumpbase \
                --split-files \
                --clip \
                --outdir TEST-HUMAN-ALT/ \
                {} ::: DRR309384 \
                       DRR309385 \
                       DRR309386 \
                       DRR309387 \
                       DRR309388 \
                       DRR309389
### Fix the automatic read mismatches between reads generated with pulling out the reads with sra-tools
echo '#!/bin/bash
f=$1
gunzip -c $f | \
sed -E "s/^((@|\+)DRR[^.]+\.[^.]+)\.(1|2)/\1/" | \
gzip -c > ${f%.fastq.gz*}-fixed.fastq.gz
' > fix_paired_end_read_names.sh
chmod +x fix_paired_end_read_names.sh
time \
parallel ./fix_paired_end_read_names.sh {} ::: $(ls TEST-HUMAN-ALT/*.fastq.gz)



echo "##################################"
echo "###                            ###"
echo "### Download reference genomes ###"
echo "###                            ###"
echo "##################################"

echo "########################################"
echo "### Download Drosophila reference genome"
echo "Release 6 plus ISO1 mitochondrial genome"
### Interactive download also availble in: https://www.ncbi.nlm.nih.gov/assembly/GCF_000001215.4#/def
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/215/GCF_000001215.4_Release_6_plus_ISO1_MT/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna.gz
gunzip GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna.gz
mv  GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna \
        Drosophila/Drosophila_reference.fasta

echo "###################################"
echo "### Download Human reference genome"
echo "Genome Reference Consortium Human Build 38 patch release 13 (GRCh38.p13)"
### Interactive download also availble in: https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.39
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.fna.gz
gunzip GCF_000001405.39_GRCh38.p13_genomic.fna.gz
mv GCF_000001405.39_GRCh38.p13_genomic.fna \
        Human/Human_reference.fasta



echo "################################################"
echo "###                                          ###"
echo "### Align Pool-seq reads into refence genome ###"
echo "###                                          ###"
echo "################################################"
### Install Heng Li's BWA and SAMTOOLS
sudo apt install -y bwa samtools
### Create index files for the reference genomes
for r in Drosophila Human
do
        cd ${r}*/
        bwa index -p ${r}_reference -a bwtsw ${r}_reference.fasta
        samtools faidx ${r}_reference.fasta
        cd -
done
### Create alignment and sorting script
echo -e '#!/bin/bash
FNAME_REF=$1
MAPQ=$2
FNAME_READ1=$3
FNAME_READ2=$4
FNAME_OUT=$(echo $FNAME_READ1 | sed s/.fastq.gz//g)
bwa mem ${FNAME_REF} ${FNAME_READ1} ${FNAME_READ2} | \
    samtools view -q ${MAPQ} -b | \
    samtools sort > ${FNAME_OUT}.bam
' > align.sh
chmod +x align.sh

echo "#####################################"
echo "### Align single-end Drosophila reads"
time \
parallel \
        ./align.sh \
                Drosophila/Drosophila_reference \
                40 \
                {1} \
                ::: $(ls Drosophila/*.fastq.gz)

echo "################################"
echo "### Align paired-end Human reads"
time \
parallel --link \
        ./align.sh \
                Human/Human_reference \
                40 \
                {1} \
                {2} \
                ::: $(ls TEST-HUMAN-ALT/*_pass_1-fixed.fastq.gz) \
                ::: $(ls TEST-HUMAN-ALT/*_pass_2-fixed.fastq.gz)


echo "#############################"
echo "###                       ###"
echo "### Pileup the alignments ###"
echo "###                       ###"
echo "#############################"
### List the bam files
for d in Drosophila Human
do
        ls ${d}/*.bam > ${d}/${d}_bam_list.txt
done
### Pileup the bam files
time \
parallel \
        samtools mpileup \
        -b {1}/{1}_bam_list.txt \
        -d 100000 \
        -q 40 \
        -Q 40 \
        -f {1}/{1}_reference.fasta \
        -o {1}/{1}.mpileup \
        ::: Drosophila Human

