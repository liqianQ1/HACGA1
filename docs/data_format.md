## Input Format

The pipeline requires a set of tab-delimited configuration files that describe the input datasets.
 All input files must be accessible inside the Docker container via the mounted data directory (default: `/data`).

### 1. Genome File (`genome.tab`)

This file specifies the reference genome sequence.

**Format:**

| Column    | Description                                    |
| --------- | ---------------------------------------------- |
| genome_id | Unique genome identifier                       |
| fasta     | Absolute or relative path to genome FASTA file |

**Example:**

```
genome1    /data/genome/genome.fa
```

------

### 2. Homolog Protein File (`homolog.tab`)

Protein homology evidence used for gene prediction.

**Format:**

| Column        | Description                |
| ------------- | -------------------------- |
| species       | Species name               |
| protein_fasta | Path to protein FASTA file |

**Example:**

```
arabidopsis    /data/protein/arabidopsis.fa
```

------

### 3. RNA-seq Evidence File (`RNAseq.tab`)

RNA-seq alignment or assembled transcript evidence.

**Format:**

| Column     | Description             |
| ---------- | ----------------------- |
| sample     | Sample name             |
| bam_or_gff | Path to BAM or GFF file |

**Example:**

```
RNAseq1    /data/rnaseq/sample1.bam
```

------

### 4. EST Evidence File (`EST.tab`)

Expressed Sequence Tag (EST) data used as additional transcription evidence.

**Format:**

| Column    | Description            |
| --------- | ---------------------- |
| est_id    | EST dataset identifier |
| est_fasta | Path to EST FASTA file |

**Example:**

```
EST1    /data/est/est.fa
```

------

## Output Details

All output files are written to the directory specified by `--Outputdir` (default: `/data`).

### 1. Gene Prediction Results

| File/Directory   | Description                                    |
| ---------------- | ---------------------------------------------- |
| `genes.gff3`     | Final integrated gene annotation (GFF3 format) |
| `genes.fasta`    | Predicted gene nucleotide sequences            |
| `proteins.fasta` | Predicted protein sequences                    |

------

### 2. Intermediate Evidence Files

| Directory     | Description                       |
| ------------- | --------------------------------- |
| `01_homolog/` | Homology-based prediction results |
| `02_rnaseq/`  | RNA-seq supported gene models     |
| `03_est/`     | EST-supported predictions         |

------

### 3. Augustus Training Models

| Path               | Description                              |
| ------------------ | ---------------------------------------- |
| `augustus_models/` | Species-specific trained Augustus models |

> ⚠️ **Important**
>  If not explicitly saved to a mounted directory, trained Augustus models will be lost when the container is removed.

### 