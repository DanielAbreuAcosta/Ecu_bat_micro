# pod5 files can give dorado the type of bascalling model to use, 
# but for fastq files this information dosent appear, so it can be good to specify the model
dorado download --model dna_r10.4.1_e8.2_260bps_fast@v4.0.0

# run basecaller and output it into a new bam file

dorado basecaller dna_r10.4.1_e8.2_260bps_fast@v4.0.0 dna_r10.4.1_e8.2_400bps_4khz-FLO_FLG114-SQK_PCB114_24-4000.pod5 > call.bam 

# alternatively, dorado can output to a fastq file using this command
dorado basecaller dna_r10.4.1_e8.2_260bps_fast@v4.0.0 dna_r10.4.1_e8.2_400bps_4khz-FLO_FLG114-SQK_PCB114_24-4000.pod5 --emit-fastq > call.fastq

# chech the first few lines of the bam file using samtools

samtools head -n 100 call.bam

#analyze all the reads for their read lengths

samtools stats call.bam | grep ^RL | cut -f 2-
  
  
#####

for i in $(seq -w 01 24)
do
samtools fastq dorado_sup_out/barcode${i}.bam > samtools_fastq_out/barcode01.fastq
   
for i in $(seq -w 01 24)
do
    barcode="barcode${i}"
    samtools fastq dorado_sup_out/barcode${i}.bam > samtools_fastq_out/barcode${i}.fastq
done


for i in $(seq -w 01 24)
  do
    NanoPlot -t 8 --fastq samtools_fastq_out/barcode${i}.fastq --loglength --plots kde --title barcode${i} -o nanoplot_out/barcode${i}
done


for i in $(seq -w 01 24)
  do
chopper --threads 8 -q 20 -l 500 -i samtools_fastq_out/barcode${i}.fastq > chopper_out/barcode${i}_filtered.fastq
done



#####

REference filter

#install ncbi installer
conda install -c conda-forge ncbi-datasets-cli

#download genome, youll need a good internet connection for this, as it will dowload a large file, 
#if it wont download you can try downloading the "dehydrated" genome and then hydrating it

datasets download genome accession GCF_022682495.2 --include genome,seq-report

# translate from bam to fastq for the host filtering (already done)
for i in $(seq -w 01 24)
  do
    barcode="barcode${i}"
    samtools fastq dorado_sup_out/barcode${i}.bam > samtools_fastq_out/barcode${i}.fastq
done

#decompress genome (might not be needed)
unzip ncbi_dataset.zip -d datasets_out

#do the host filtering
for i in $(seq -w 01 24)
  do
    barcode="barcode${i}"
    minimap2 -a -x lr:hq datasets_out/ncbi_dataset/data/GCF_022682495.2/GCF_022682495.2_HLdesRot8A.1_genomic.fna samtools_fastq_out/barcode${i}.fastq | samtools view -bh | samtools sort -o minimap_out/barcode${i}.bam
done

#filter the non-host sequences
for i in $(seq -w 01 24)
  do
    barcode="barcode${i}"
samtools fastq -f 4 minimap_out/barcode${i}.bam | samtools sort -o samtools_filter_out/barcode${i}.fastq
done

    minimap2 -a -x lr:hq datasets_out/ncbi_dataset/data/GCF_022682495.2/GCF_022682495.2_HLdesRot8A.1_genomic.fna samtools_fastq_out/barcode04.fastq | samtools view -bh | samtools sort -o minimap_out/barcode04.bam
samtools fastq -f 4 minimap_out/barcode04.bam | samtools sort -o samtools_filter_out/barcode04.fastq


If you have to stop a command and continue it later, you can hit Ctl+z to pause it, and input fg to contiue when you can

