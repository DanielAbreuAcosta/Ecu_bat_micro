# Bioinformatics Pipeline for UTI Pathogen

**DISCLAIMER: this code is nowhere near finished**

*Identification from Host-Depleted Nanopore Sequencing Data (Optimized for bacterial detection in bat urinary samples/tissues relevant to UTI-causing pathogens)*

*Using an Apple Silicon Macbook (M2 chip), some of this code will still work in other operating systems, especially the Python code. Most, if not all, of the installation steps won't work on other OS*

## Step 0: Installations

### Dorado

1.  To download the latest version of Dorado, go to <https://nanoporetech.com/software/other/dorado> and select your operating system. If you don't have an account with ONT, you'll be prompted to create before installing.

2.  Decompress the downloaded zip file

3.  Open terminal and navigate to the directory containing the file

4.  To access dorado from anywhere on the computer, move it to the local directory and add it to the configuration file:

``` {.console style="color: gray"}
sudo mv dorado-0.9.1-osx-arm64 /usr/local/dorado

echo 'export PATH=/usr/local/dorado/bin:$PATH' >> ~/.zshrc

source ~/.zshrc
```

5.  run `dorado --version` to check that the installation was successful, this should work from any directory

### Homebrew

### Samtools

### Conda

### Bioconda

### NanoPlot

## Step 1: Basecalling & Quality Control

### Basecalling

1.  Loop through each barcode directory and convert raw Nanopore signals into DNA sequences using Dorado:

``` {.console style="color: gray"}
for i in $(seq -w 01 24)  # -w ensures leading zeros
do
    barcode="barcode${i}"
    dorado basecaller dna_r10.4.1_e8.2_400bps_sup@v5.0.0 "pod5_pass/$barcode/" > "dorado_sup_out/${barcode}.bam"
done
```

`dna_r10.4.1_e8.2_400bps_sup@v5.0.0` is the version of the basecaller to use, pod5 can give dorado the kind of basecalling model to use, but it can be good to specify the model, specially if using fastq files (which don't contain this information).

2.  To check any of the BAM files, run:

``` {.console style="color: gray"}
samtools head -n 100 barcodexx.bam
```

Remember to replace the `barcodexx.bam` with your filename.

### Data Visualization

3.  In order to clean up the data, first look into the quality of the reads, as well as their length, this can be done with NanoPlot.

This package prefers fastq files, so first transform them using samtools:

``` {.console style="color: gray"}
samtools fastq dorado_sup_out/barcode01.bam > samtools_fastq_out/barcode01.fastq
```

4.  Now, use NanoPlot to generate a report for each barcode:

``` {.conda style="color: green"}
NanoPlot -t 2 --fastq samtools_fastq_out/barcode01.fastq --loglength --plots kde --title barcode01 -o nanoplot_out/barcode01
```

`-t 2` makes the package run in two threads, this can speed up the process with more powerful gpus.

`--fastq samtools_fastq_out/barcode01.fastq` gives the file type and the file name.

`--loglength` puts read lengths on a logarithmic scale.

`kde` creates a kernel density estimate plot (kde plot), which is a method for visualizing the distribution of the data, similar to a histogram (<https://seaborn.pydata.org/generated/seaborn.kdeplot.html>).

`--title barcode01` Puts the title on all plots.

`-o nanoplot_out/barcode01` puts all generated graphs in this directory.

\_\_\_\_

### Quality Filtering

Remove short and low-quality reads using Filtlong or NanoFilt:

``` console
bash
Copy
Edit
filtlong --min_length 1000 --keep_percent 90 basecalled.fastq > filtered.fastq
```

## Step 2: Taxonomic Classification (UTI Pathogen Focus)

### Bacterial Identification (Including UTI Pathogens)

Use Kraken2, Centrifuge, or Kaiju for taxonomic assignment:

``` console
bash
Copy
Edit
kraken2 --db kraken_db --threads 8 --report bacteria_report.txt --output bacteria_classified.txt --use-names filtered.fastq
```

#### Key UTI Pathogens to Check For:

-   *Escherichia coli* (UPEC - Uropathogenic *E. coli*)
-   *Klebsiella pneumoniae*
-   *Proteus mirabilis*
-   *Enterococcus faecalis*
-   *Staphylococcus saprophyticus*
-   *Pseudomonas aeruginosa*
-   *Morganella morganii*

### Visualization of UTI Pathogens

Generate interactive taxonomic plots using KronaTools.

``` console
bash
Copy
Edit
cut -f2,3 bacteria_classified.txt | ktImportTaxonomy -o bacteria_krona.html
```

## Step 3: UTI Pathogen Confirmation & Functional Annotation

### Targeted UTI Pathogen Screening

Use MetaPhlAn or PathoScope to refine UTI pathogen detection.

``` console
bash
Copy
Edit
metaphlan filtered.fastq --input_type fastq -o pathogen_abundance.txt
```

Alternative: Check for pathogenic UTI genes using MASH screen.

``` console
bash
Copy
Edit
mash screen -w -i 0.9 UTI_reference_db.msh filtered.fastq > mash_results.txt
```

### Antimicrobial Resistance (AMR) Screening for UTI Bacteria

Use ResFinder or CARD (Comprehensive Antibiotic Resistance Database).

``` console
bash
Copy
Edit
rgi main --input_sequence filtered.fastq --output rgi_output.txt --aligner DIAMOND --local
```

#### Look for AMR genes common in UTI pathogens, such as:

-   **Beta-lactam resistance**: *blaCTX-M, blaTEM, blaSHV*
-   **Fluoroquinolone resistance**: *gyrA, parC*
-   **Aminoglycoside resistance**: *aac(6')-Ib*
-   **Sulfonamide resistance**: *sul1, sul2*

### Virulence Factor Detection for UTI Pathogens

Use ABRICATE with the VFDB (Virulence Factor Database).

``` console
bash
Copy
Edit
abricate --db vfdb filtered.fastq > virulence_report.txt
```

#### Key virulence genes in UTI bacteria to look for:

-   *Escherichia coli* (*fimH, papG, sfa, iroN*)
-   *Proteus mirabilis* (*hpmA, mrpA*)
-   *Klebsiella pneumoniae* (*rmpA, yersiniabactin*)

## Step 4: Bacterial Genome Assembly & UTI Pathogen Strain Analysis

### De Novo Assembly of UTI Pathogen Genomes

Use Flye for assembling long-read bacterial genomes.

``` console
bash
Copy
Edit
flye --nano-raw filtered.fastq --out-dir assembly_output --genome-size 5m
```

### Polishing Assembly to Improve Accuracy

Use Medaka for error correction.

``` console
bash
Copy
Edit
medaka_consensus -i filtered.fastq -d assembly_output/assembly.fasta -o polished_assembly
```

### Confirm UTI Pathogen Identity via BLAST

Align assembled genomes to NCBI’s bacterial reference database.

``` console
bash
Copy
Edit
blastn -query polished_assembly/consensus.fasta -db nt -out blast_results.txt -outfmt 6
```

## Step 5: Phylogenetics & Comparative Genomics for UTI Pathogens

### Phylogenetic Placement of UTI Bacteria

Use Mashtree to analyze evolutionary relationships.

``` console
bash
Copy
Edit
mashtree --numcpus 4 polished_assembly/*.fasta > phylogeny_tree.nwk
```

### Comparative Genomics for UTI Pathogens

Use Roary or Panaroo for pan-genome analysis.

``` console
bash
Copy
Edit
roary -e -n -v *.gff
```

## Key Additions for UTI Pathogen Analysis

✅ Focus on key UTI bacteria (*E. coli, Klebsiella, Proteus*, etc.)\
✅ Detect antimicrobial resistance genes relevant to UTI treatment\
✅ Identify virulence genes involved in UTI pathogenesis\
✅ Ensure accurate pathogen strain identification

## Final Outputs

-   List of UTI pathogens present in the sample\
-   AMR profile of identified bacteria\
-   Virulence factor annotations\
-   Assembled bacterial genomes\
-   Phylogenetic relationships of UTI-causing bacteria
