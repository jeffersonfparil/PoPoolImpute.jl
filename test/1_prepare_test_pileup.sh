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
mkdir Drosophila-Wei-etal-2017/ ### https://doi.org/10.1534/genetics.116.197335
### Note: fasterq-dump is the new software set to replace fastq-dump, but in my opinion it is still immature as it lacks some of the features that fastq-dump has, e.g. --gzip
time parallel fastq-dump \
                --gzip \
                --skip-technical \
                --readids \
                --read-filter pass \
                --dumpbase \
                --split-files \
                --clip \
                --outdir Drosophila-Wei-etal-2017/ \
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
### Note: 100 bp paired-end reads
mkdir Human-DiNatale-etal-2021/ ### https://dx.doi.org/10.52733%2FKCJ18n2.a1
time parallel fastq-dump \
                --gzip \
                --skip-technical \
                --readids \
                --read-filter pass \
                --dumpbase \
                --split-files \
                --clip \
                --outdir Human-DiNatale-etal-2021/ \
                {} :::  SRR11801595 SRR11801596 SRR11801597 SRR11801598 SRR11801599 \
                        SRR11801600 SRR11801601 SRR11801602 SRR11801603 SRR11801604 \
                        SRR11801605 SRR11801606 SRR11801607 SRR11801608 SRR11801609 \
                        SRR11801610 SRR11801611 SRR11801612 SRR11801613 SRR11801614 \
                        SRR11801615 SRR11801616 SRR11801617 SRR11801618 SRR11801619 \
                        SRR11801620 SRR11801621 SRR11801622 SRR11801623 SRR11801624 \
                        SRR11801625 SRR11801626 SRR11801627 SRR11801628 SRR11801629 \
                        SRR11801630 SRR11801631 SRR11801632 SRR11801633 SRR11801634 \
                        SRR11801635 SRR11801636 SRR11801637 SRR11801638 SRR11801639 \
                        SRR11801640 SRR11801641 SRR11801642



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
        Drosophila-Wei-etal-2017/Drosophila_reference.fasta

echo "###################################"
echo "### Download Human reference genome"
echo "Genome Reference Consortium Human Build 38 patch release 13 (GRCh38.p13)"
### Interactive download also availble in: https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.39
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.fna.gz
gunzip GCF_000001405.39_GRCh38.p13_genomic.fna.gz
mv GCF_000001405.39_GRCh38.p13_genomic.fna \
        Human-DiNatale-etal-2021/Human_reference.fasta



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
                Drosophila-Wei-etal-2017/Drosophila_reference \
                40 \
                {1} \
                ::: $(ls Drosophila-Wei-etal-2017/*.fastq.gz)

echo "################################"
echo "### Align paired-end Human reads"
time \
parallel --link \
        ./align.sh \
                Human-DiNatale-etal-2021/Human_reference \
                40 \
                {1} \
                {2} \
                ::: $(ls Human-DiNatale-etal-2021/*.fastq.gz | grep _pass_1) \
                ::: $(ls Human-DiNatale-etal-2021/*.fastq.gz | grep _pass_2)



echo "#############################"
echo "###                       ###"
echo "### Pileup the alignments ###"
echo "###                       ###"
echo "#############################"
### List the bam files
for DIR in Drosophila-Wei-etal-2017 Human-DiNatale-etal-2021
do
        ls ${DIR}/*.bam > ${DIR}/Drosophila_bam_list.txt
done
### Pileup the bam files
time \
parallel \
        samtools mpileup \
        -b {1}*/{1}_bam_list.txt \
        -d 100000 \
        -q 40 \
        -Q 40 \
        -f {1}_reference \
        -o {1}*/{1}.mpileup \
        ::: Drosophila Human
