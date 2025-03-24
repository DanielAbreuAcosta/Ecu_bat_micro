# Bioinformatics Pipeline for UTI Pathogen

*Identification from Host-Depleted Nanopore Sequencing Data (Optimized for bacterial detection in bat urinary samples/tissues relevant to UTI-causing pathogens)*

## Step 1: Basecalling & Quality Control

### Basecalling

Convert raw Nanopore signals into DNA sequences using Guppy or Bonito.

``` console

# pod5 files can give dorado the type of bascalling model to use, 
# but for fastq files this information dosent appear, so it can be good to specify the model

dorado download --model dna_r10.4.1_e8.2_260bps_fast@v4.0.0

# run basecaller and output it into a new bam file

dorado basecaller dna_r10.4.1_e8.2_260bps_fast@v4.0.0 dna_r10.4.1_e8.2_400bps_4khz-FLO_FLG114-SQK_PCB114_24-4000.pod5 > call.bam

# chech the first few lines of the bam file using samtools

samtools head -n 100 call.bam
```

### Quality Filtering

Remove short and low-quality reads using Filtlong or NanoFilt.

``` console
bash
Copy
Edit
filtlong --min_length 1000 --keep_percent 90 basecalled.fastq > filtered.fastq
```

## Step 2: Taxonomic Classification (UTI Pathogen Focus)

### Bacterial Identification (Including UTI Pathogens)

Use Kraken2, Centrifuge, or Kaiju for taxonomic assignment.

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
