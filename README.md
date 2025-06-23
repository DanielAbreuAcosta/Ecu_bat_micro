# Bioinformatics Pipeline for Bat Pathogens

**DISCLAIMER: this code is nowhere near finished**

*Identification from Host-Depleted Nanopore Sequencing Data (Optimized for bacterial detection of known human pathogens in bat tissues).*

*Done using an Apple Silicon Macbook (M2 chip), some of this code will still work in other operating systems, especially the Python code. Most, if not all, of the installation steps won't work on other OS.*

## Step 0: Installations

### Dorado

A basecalling algorithm developed by ONT for decoding the raw electrical signals produced from sequencing, and translating them into nucleotide sequences. It also outputs a quality score for each base, which will be used later.

1.  To download the latest version of Dorado, go to <https://nanoporetech.com/software/other/dorado> and select your operating system. If you don't have an account with ONT, you'll be prompted to create one before installing.

2.  Decompress the downloaded zip file

3.  Open terminal and navigate to the directory containing the file

4.  To access Dorado from anywhere on the computer, move it to the local directory and add it to the configuration file:

``` bash
sudo mv dorado-0.9.1-osx-arm64 /usr/local/dorado

echo 'export PATH=/usr/local/dorado/bin:$PATH' >> ~/.zshrc

source ~/.zshrc
```

5.  run `dorado --version` to check that the installation was successful, this should work from any directory

### Homebrew

An open-source package manager for macOS and Linux, it will help with the installation of packages like Samtools.

1.  To install Homebrew, run:

``` bash
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```

2.  to check that installation was sucesfull, run `brew config`.

### Samtools

A program for interacting with high-throughput sequencing data, it is useful for reading, writing, editing, indexing and viewing data,

1.  **After installing homebrew**, run:

``` bash
brew install samtools
```

2.  Check that installation was successful by running `samtools version`.

### Conda

Anaconda is a program that aids in Python programming, allowing for easy download of packages, it also includes Conda, which is an open-source manager.

1.  To install it, go to <https://www.anaconda.com/download/success>, and select your operating system.
2.  Execute the downloaded file, and follow the instructions of the installer.
3.  Once downloaded, go to terminal. If installed correctly you should see `(base)` before the user and location information, this indicates that you are now in the base environment of conda, instead of the normal terminal
4.  To exit out of conda run `conda deactivate`, and the `(base)` tag should dissapear. You are now back to the normal terminal. To go back into conda run `conda activate`.

#### Set up an environment

`(base)` indicates that you are in the base environment of conda, this contains all available packages. Best practice is to never work in this base environment, but to make a new environment for each new project instead.

1.  To create an environment, run `conda create --name example`.
2.  Check all current environments by running: `conda env list`.
3.  Activate an environment with: `conda activate example`.
4.  See all packages inside the environment with: `conda list` (some packages come preinstalled).
5.  To exit back to the base environment, run: `conda deactivate`.

### Bioconda

Bioconda is a package manager (or channel) that lets you install software packages related to biomedical research.

1.  In the base environment, download Bioconda, as well as Conda Forge (needed for Bioconda to work)

``` python
conda config --add channels bioconda
conda config --add channels conda-forge
```

2.  Make Conda strictly follow the order of channels when resolving dependencies, meaning it will first use the packages from bioconda and then the ones from conda forge. This prevents conda from looking for the same package in multiple channels.

``` python
conda config --set channel_priority strict
```

### Minimap

``` python
conda activate example
conda install -c bioconda minimap2
```

### 

### NanoPlot

Nanoplot is a plotting tool for long read sequencing data and alignments. Go to your working environment and install it:

``` python
conda activate example
conda install -c bioconda nanoplot
```

### Chopper

Chopper is a tool to filter and trim fastq files from long read sequencing technologies like ONT. Go to your working environment and install it:

``` python
conda activate example
conda install -c bioconda chopper
```

### Flye

Flye is a contig assembler developed for long read sequencing. We will mainly be using the metaFlye algorithm withing the Flye library, which is intended assembling metagenomic data with highly uneven coverage (Milkhail et al., 2020) (<https://doi.org/10.1038/s41592-020-00971-x>), and has been shown to work better that other assemblers (Latorre‑Pérez et al., 2020) (<https://doi.org/10.1038/s41598-020-70491-3>). To install Flye:

``` python
conda activate example
conda install -c bioconda flye
```

## Step 1: Basecalling & Quality Control

### Basecalling

1.  Loop through each barcode directory and convert raw Nanopore signals into DNA sequences using Dorado (WARNING: this step is very taxing on computer resources, it will take a long time to run on a normal computer. It is probably best to do this step on a cluster if you have access to one):

``` bash
for i in $(seq -w 01 24)  # -w ensures leading zeros
do
    barcode="barcode${i}"
    dorado basecaller dna_r10.4.1_e8.2_400bps_sup@v5.0.0 "pod5_pass/$barcode/" > "dorado_sup_out/${barcode}.bam"
done
```

`dna_r10.4.1_e8.2_400bps_sup@v5.0.0` is the version of the basecaller to use, pod5 can give dorado the kind of basecalling model to use, but it can be good to specify the model, specially if using fastq files (which don't contain this information).

Each individual base-call is accompanied by another character which indicate the error probability for that base-call, called the Phred Quality score, which is calculated with the following formula:

Q = −10 × log10 p

Where p represents the estimated error probability. For example, a base-call with a Q score of 20 would have a probability of 1/100 of being incorrect. For a more in-depth explanation on this, read Ewing & Green (1998) (<http://genome.cshlp.org/cgi/pmidlookup?view=long&pmid=9521922>)

2.  To check any of the generated BAM files, run:

``` bash
samtools head -n 100 barcodexx.bam
```

Remember to replace the `barcodexx.bam` with your file-name.

add the bash script to run in the cluster. Add notes

### Host Sequence Filtering

Even though our data was generated using host-depleted sequencing, the host genome can still make it through. Because of this, it is necessary to remove these sequences before proceeding. For this we will be using the latest reference genome of *Desmodus rotundus*.

leave this as an example for one and then add the real process for multiple host species

1.  download the reference genome from NCBI

``` bash
datasets download genome
accession GCF_022682495.2 --include
gff3,rna,cds,protein,genome,seq-report
```

2.  filter using minimap and samtools <https://linsalrob.github.io/ComputationalGenomicsManual/Deconseq/>

``` bash
minimap2 --split-prefix=tmp$$ -a -xsr GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz R1.fastq.gz R2.fastq.gz | samtools view -bh | samtools sort -o output.bam
```

### Data Visualization

3.  In order to clean up the data, first look into the quality of the reads, as well as their length, this can be done with NanoPlot. This package prefers fastq files, so first transform them using samtools:

``` bash
for i in $(seq -w 01 24)
  do
    barcode="barcode${i}"
    samtools fastq dorado_sup_out/barcode${i}.bam > samtools_fastq_out/barcode${i}.fastq
done
```

4.  Now, use NanoPlot to generate a report for each barcode:

``` python
for i in $(seq -w 01 24)
  do
    NanoPlot -t 8 --fastq samtools_fastq_out/barcode${i}.fastq --loglength --plots kde --title barcode${i} -o nanoplot_out/barcode${i}
done
```

`-t 8` makes the package run in eight threads, this can speed up the process with more powerful gpus.

`--fastq samtools_fastq_out/barcode${i}.fastq` gives the file type and the file name.

`--loglength` puts read lengths on a logarithmic scale.

`kde` creates a kernel density estimate plot (kde plot), which is a method for visualizing the distribution of the data, similar to a histogram (<https://seaborn.pydata.org/generated/seaborn.kdeplot.html>).

`--title barcode${i}` Puts the title on all plots.

`-o nanoplot_out/barcode${i}` puts all generated graphs in this directory.

### Data Cleaning

5.  Remove short and low-quality reads using chopper:

``` python
for i in $(seq -w 01 24)
  do
chopper --threads 8 -q 20 -l 20 -i samtools_filter_out/barcode${i}.fastq > chopper_out/barcode${i}_filtered.fastq
done
```

`--threads 8` makes the package run in eight threads, this can speed up the process with more powerful gpus.

`-q 20` sets the minimum Phred Quality Score to 20. Any sequence with an average score below this will get removed.

`-l 20` sets the minimum read length to 20bp.

`-i samtools_fastq_out/barcode${i}.fastq` tells the package the name of the file to filter.

## Step 2: Taxonomic Classification (UTI Pathogen Focus)

### De-novo Assembly

1.  Assemble contigs using metaFlye (WARNING: this step is relatively taxing on computer resources, it might take a while to run the command):

``` python
for i in $(seq -w 01 24)
  do
    barcode="barcode${i}"
flye --meta --read-error 0.03 --nano-hq chopper_out/barcode${i}_filtered.fastq --out-dir flye_out/barcode${i} --threads 8
done
```

`--meta` activates metaFlye mode.

`--read-error 0.03` specifies that data has already been filtered to Q20 (in step 1.5.).

`--nano-hq` specifies that the data was basecalled using Dorado's super accurate mode (in step 1.1.).

`chopper_out/barcode${i}_filtered.fastq` tells the package the name of the file to assemble.

`--out-dirflye_out/barcode${i}` specifies the output directory

`--threads 8` makes the package run in eight threads, this can speed up the process with more powerful gpus.

2.  Extract the assemblies from their enclosing folders:

``` python
for i in $(seq -w 01 24)
  do
    barcode="barcode${i}"
    assembly="assembly{i}"
cp flye_out/barcode${i}/assembly.fasta flye_assembly/assembly${i}.fasta
cp flye_out/barcode${i}/assembly_info.txt flye_assembly/assembly${i}_info.txt
done
```

### UTI Identification

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
