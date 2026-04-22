#!/bin/bash

cd TrimGalore-0.6.7
chmod +x trim_galore
export PATH=$(pwd):$PATH
cd ..
trim_galore --version
chmod +x FastQC/fastqc
export PATH=$(pwd)/FastQC:$PATH
fastqc --version
cd hisat2-2.2.1
echo 'export PATH=$PWD:$PATH' >> ~/.bashrc
source ~/.bashrc
cd ..
hisat2 --version

SPECIES=plutella
SPECIES2=plutella
FULL=Plutella_xylostella

echo "${FULL}"

mkdir -p ${SPECIES}_fastq ${SPECIES}_trimmed ${SPECIES}_mapped logs


### Sanity checks – required files
set -euo pipefail

echo "Running sanity checks..."

# SRR list
if [ ! -f ${SPECIES}_SRR.txt ]; then
    echo "ERROR: ${SPECIES}_SRR.txt not found"
    exit 1
fi

# Annotation
if [ ! -f ${FULL}.gff3 ]; then
    echo "ERROR: Annotation file ${FULL}.gff3 not found"
    exit 1
fi

# HISAT2 index
if [ ! -f ${SPECIES2}_index_hisat/${SPECIES2}_index.1.ht2 ]; then
    echo "ERROR: HISAT2 index not found in ${SPECIES2}_index_hisat/"
    echo "Expected files like: ${SPECIES2}_index.1.ht2"
    exit 1
fi

# Required tools
for tool in fasterq-dump trim_galore fastqc hisat2 samtools featureCounts gffread; do
    if ! command -v $tool >/dev/null 2>&1; then
        echo "ERROR: $tool not found in PATH"
        exit 1
    fi
done

echo "All sanity checks passed"
echo

while read -r SRR; do
    [[ -z "$SRR" ]] && continue
    echo "Processing $SRR"

    #1. Downloading
    echo "Downloading $SRR"
    DOWNLOAD_OK=0
    for i in {1..5}; do
        echo "Attempt $i for $SRR"
        fasterq-dump \
          -O ${SPECIES}_fastq \
          --progress \
          "$SRR" && DOWNLOAD_OK=1 && break
        sleep 30
    done

    if [ $? -ne 0 ]; then
        echo "ERROR: Download failed for $SRR"
        echo "Skipping $SRR"
        continue
    fi

    #2. Trimming
    file1=${SPECIES}_fastq/${SRR}_1.fastq
    file2=${SPECIES}_fastq/${SRR}_2.fastq

    if [[ -f "$file1" && -f "$file2" ]]; then
        echo "trimming paired-end reads for $SRR"
        TrimGalore-0.6.7/trim_galore --paired \
            --quality 20 \
            --fastqc \
            --length 36 \
            -o ${SPECIES}_trimmed \
            "$file1" "$file2" 2>&1 | tee logs/${SRR}_trim.log
        if [ $? -eq 0 ]; then
            rm "$file1" "$file2"
        fi
    else
        echo "Trimming single-end reads for $SRR"
        TrimGalore-0.6.7/trim_galore \
            --quality 20 \
            --fastqc \
            --length 36 \
            -o ${SPECIES}_trimmed \
            "${SPECIES}_fastq/${SRR}.fastq" 2>&1 | tee logs/${SRR}_trim.log
        if [ $? -eq 0 ]; then
            rm "${SPECIES}_fastq/${SRR}.fastq"
        fi
    fi

    #3. Mapping
    #paired-end mapping
    if ls ${SPECIES}_trimmed/${SRR}_1_val_1.fq 1>/dev/null 2>&1; then
        echo "Mapping paired-end reads for $SRR"
        #--large-index for locusta with .ht2l files
        hisat2 --dta \
            -x ${SPECIES2}_index_hisat/${SPECIES2}_index \
            -1 ${SPECIES}_trimmed/${SRR}_1_val_1.fq \
            -2 ${SPECIES}_trimmed/${SRR}_2_val_2.fq \
            -S ${SPECIES}_mapped/${SRR}.sam 2>&1 | tee logs/${SRR}_map.log

        samtools view -bS ${SPECIES}_mapped/${SRR}.sam | samtools sort -o ${SPECIES}_mapped/${SRR}.bam
        rm ${SPECIES}_mapped/${SRR}.sam ${SPECIES}_trimmed/${SRR}_1_val_1.fq ${SPECIES}_trimmed/${SRR}_2_val_2.fq
    fi

    #single-end mapping
    if ls ${SPECIES}_trimmed/${SRR}_trimmed.fq 1>/dev/null 2>&1; then
        echo "Mapping single-end reads for $SRR"
        hisat2 --dta \
            -x ${SPECIES2}_index_hisat/${SPECIES2}_index \
            -U ${SPECIES}_trimmed/${SRR}_trimmed.fq \
            -S ${SPECIES}_mapped/${SRR}.sam 2>&1 | tee logs/${SRR}_map.log

        samtools view -bS ${SPECIES}_mapped/${SRR}.sam | samtools sort -o ${SPECIES}_mapped/${SRR}.bam
        rm ${SPECIES}_mapped/${SRR}.sam ${SPECIES}_trimmed/${SRR}_trimmed.fq
    fi

    echo "Finished $SRR"
done < ${SPECIES}_SRR.txt

#Convert GFF3 to GTF
if [ ! -f ${FULL}.gtf ]; then
    echo "Converting GFF3 to GTF"
    gffread ${FULL}.gff3 -T -o ${FULL}.gtf
fi

#4. Feature counting BAM
#-p for paired reads
echo "Counting features"
export PATH=$PWD/subread-2.0.0-Linux-x86_64/bin:$PATH
featureCounts -p \
  -t exon \
  -g gene_id \
  -O \
  -M \
  -a ${FULL}.gtf \
  -o ${SPECIES}_counts.txt \
  ${SPECIES}_mapped/*.bam 2>&1 | tee logs/featureCounts_${SPECIES}.log

echo "All finished"
        

