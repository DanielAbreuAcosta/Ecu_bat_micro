#!/bin/bash

# before running make sure that the following files are present: mapping.tsv, all the barcodeXX.fastq files contained in samtools_fastq_out

#create working directories
mkdir -p datasets_out/
mkdir -p minimap_out/
mkdir -p samtools_filter_out/
mkdir -p samtools_fastq_out/

#create working environment
conda create --name ecu_bat_micro
conda activate ecu_bat_micro

#download bioconda and conda forge
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict

#install ncbi installer, minimap2 and samtools
conda install -c conda-forge ncbi-datasets-cli -c bioconda minimap2 samtools

#download reference genomes
datasets download genome accession GCA_004027475.1 --include genome,seq-report --filename GCA_004027475.1.zip
unzip GCA_004027475.1.zip -d datasets_out/GCA_004027475.1

datasets download genome accession GCA_004027735.1 --include genome,seq-report --filename GCA_004027735.1.zip
unzip GCA_004027735.1.zip -d datasets_out/GCA_004027735.1

datasets download genome accession GCA_014824575.3 --include genome,seq-report --filename GCA_014824575.3.zip
unzip GCA_014824575.3.zip -d datasets_out/GCA_014824575.3

datasets download genome accession GCA_027563665.1 --include genome,seq-report --filename GCA_027563665.1.zip
unzip GCA_027563665.1.zip -d datasets_out/GCA_027563665.1

datasets download genome accession GCA_027574615.1 --include genome,seq-report --filename GCA_027574615.1.zip
unzip GCA_027574615.1.zip -d datasets_out/GCA_027574615.1

datasets download genome accession GCA_036850765.1 --include genome,seq-report --filename GCA_036850765.1.zip
unzip GCA_036850765.1.zip -d datasets_out/GCA_036850765.1

datasets download genome accession GCA_038363175.3 --include genome,seq-report --filename GCA_038363175.3.zip
unzip GCA_038363175.3.zip -d datasets_out/GCA_038363175.3

datasets download genome accession GCA_039880945.1 --include genome,seq-report --filename GCA_039880945.1.zip
unzip GCA_039880945.1.zip -d datasets_out/GCA_039880945.1

datasets download genome accession GCA_963259705.2 --include genome,seq-report --filename GCA_963259705.2.zip
unzip GCA_963259705.2.zip -d datasets_out/GCA_963259705.2

#Host Sequence Filtering
# Declare an associative array
declare -A barcode_map

# Read mapping from TSV file
while IFS=$'\t' read -r barcode _ _ _ genome_code _; do
    barcode_map[$barcode]=$genome_code
done < mapping.tsv 

# Loop through barcode numbers
for i in $(seq -w 01 24); do
    barcode="barcode${i}"
    genome_code=${barcode_map[$barcode]}
    
    # Check if the barcode has a genome code mapping
    if [ -n "$genome_code" ]; then
        echo "Processing $barcode with $genome_code"
        minimap2 -a -x lr:hq datasets_out/${genome_code}/ncbi_dataset/data/${genome_code}/${genome_code}*.fna \
            samtools_fastq_out/${barcode}.fastq | \
            samtools view -bh | \
            samtools sort -o minimap_out/${barcode}.bam
    else
        echo "No genome code found for $barcode, skipping..."
    fi
done

#filter the non-host sequences
for i in $(seq -w 01 24)
  do
    barcode="barcode${i}"
samtools fastq -f 4 minimap_out/barcode${i}.bam | samtools sort -o samtools_filter_out/barcode${i}.fastq
done