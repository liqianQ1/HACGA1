## Input Format

The pipeline requires a set of tab-delimited configuration files that describe the input datasets.

 All input files must be accessible inside the Docker container via the mounted data directory (default: `/data`).



### 1. Genome File (`genome.tab`)

This file specifies the reference genome sequence.

| Sample      | Species      | Genome_fasta           | Masked_genome_fasta                  |
| ----------- | ------------ | ---------------------- | ------------------------------------ |
| sample_name | species_name | /data/genome/genome.fa | /data/genome/genome.repeat.masked.fa |

------

### 2. Homolog Protein File (`homolog.tab`)

Protein homology evidence used for gene prediction.

This file can contain a single row specifying a shared protein database for all samples.

| Homolog_set_name | Homolog_protein_fasta      |
| ---------------- | -------------------------- |
| Homolog_name     | /data/genome/Homolog.fasta |

------

### 3. RNA-seq Evidence File (`RNAseq.tab`)

RNA-seq alignment or assembled transcript evidence.

| Sample       | Fastq1                                          | Fastq2                                          |
| ------------ | ----------------------------------------------- | ----------------------------------------------- |
| sample_name1 | /data/genome/RNA-seq/sample_name1_1.clean.fq.gz | /data/genome/RNA-seq/sample_name1_2.clean.fq.gz |
| sample_name2 | /data/genome/RNA-seq/sample_name2_1.clean.fq.gz | /data/genome/RNA-seq/sample_name2_2.clean.fq.gz |

### 4. EST Evidence File (`EST.tab`)

Expressed Sequence Tag (EST) data used as additional transcription evidence.

This file can contain a single row specifying a shared Trinity FASTA file for all samples.

| sample      | Trinity_fasta              |
| ----------- | -------------------------- |
| sample_name | /data/genome/Trinity.fasta |

------

## Output Details

All output files are written to the directory specified by `--Outputdir` .

### 1. Genome-level *Ab Initio* Prediction Results

| File                          | Description                                                  |
| ----------------------------- | ------------------------------------------------------------ |
| `03augustus/augustus.gff`     | Gene structures predicted by **Augustus** from the genome.   |
| `04genemarket/genemarket.gff` | Gene structures predicted by **GeneMark-ET** from the genome sequence. |
| `06glimmerhmm/glimmerhmm.gff` | Gene structures predicted by **GlimmerHMM** from the genome sequence. |

------

### 2. Protein-level Prediction Results

| File                          | Description                                                  |
| ----------------------------- | ------------------------------------------------------------ |
| `05genemarkep/genemarkep.gff` | Gene structures predicted by **GeneMark-EP+** from the genome. |

------

### 3. Transcript-level Prediction Results

| File                                     | Description                                                  |
| ---------------------------------------- | ------------------------------------------------------------ |
| `01RNAdenovo/03_stringtie/RNADenovo.gff` | Transcripts assembled by **StringTie** from RNA-seq data.    |
| `02trinity_pasa/pasa.1.end.gff`          | Terminal exon structures predicted by **PASA** for Trinity transcripts. |

### 4. Integrated Annotation

| File                       | Description                                                  |
| -------------------------- | ------------------------------------------------------------ |
| `07evm/evm.gff3`           | Final integrated gene models produced by **EvidenceModeler (EVM)**, combining predictions from Augustus, GeneMark-EP+,GeneMark-ET, GlimmerHMM, PASA and other evidence sources (GFF3 format). Includes both protein‑coding and non‑coding genes. |
| `08evm_pasa/*.protein.gff` | Protein‑coding gene structures refined by PASA based on EVM results, one file per processed assembly or configuration (GFF3 format). Represents high‑confidence coding gene models supported by RNA‑seq evidence. |

