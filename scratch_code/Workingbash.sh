# pod5 files can give dorado the type of bascalling model to use, 
# but for fastq files this information dosent appear, so it can be good to specify the model
dorado download --model dna_r10.4.1_e8.2_260bps_fast@v4.0.0

# run basecaller and output it into a new bam file

dorado basecaller dna_r10.4.1_e8.2_260bps_fast@v4.0.0 dna_r10.4.1_e8.2_400bps_4khz-FLO_FLG114-SQK_PCB114_24-4000.pod5 > call.bam 

# alternatively, dorado can output to a fastq file using this command
dorado basecaller dna_r10.4.1_e8.2_260bps_fast@v4.0.0 dna_r10.4.1_e8.2_400bps_4khz-FLO_FLG114-SQK_PCB114_24-4000.pod5 --emit-fastq


# chech the first few lines of the bam file using samtools

samtools head -n 100 call.bam

#analyze all the reads for their read lengths

samtools stats call.bam | grep ^RL | cut -f 2-
  


  