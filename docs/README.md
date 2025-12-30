# 🧬 HACGA1  
**An Integrated, Evidence-Driven Genome Annotation Pipeline**

HACGA1 is a fully automated genome structural annotation pipeline designed for **eukaryotic genomes**.  
It integrates *ab initio gene prediction*, *transcriptome-based evidence*, and *homology-based inference* into a unified, reproducible workflow.

The pipeline is implemented in Python and deployed exclusively via Docker, ensuring cross-platform consistency, reproducibility, and minimal environment configuration.

---

## ✨ Key Features

- End-to-end automated genome annotation
- Integration of multiple annotation strategies:
  - **Ab initio prediction**: AUGUSTUS, GlimmerHMM, GeneMark
  - **Transcript evidence**: RNA-seq, Trinity, PASA, TransDecoder
  - **Homology evidence**: GeneMark-EP+ with protein alignments
- Evidence integration and consensus modeling using **EvidenceModeler (EVM)**
- UTR and alternative splicing refinement with **PASA**
- Docker-based deployment for reproducibility
- Parallel execution with multi-threading and optional SLURM support
- Robust error handling and modular pipeline design

---

## 🧠 Pipeline Overview

HACGA1 performs genome annotation through the following major stages:

1. Genome preprocessing and sequence partitioning for parallel computation
2. RNA-seq alignment and transcript assembly (HISAT2 + StringTie)
3. De novo transcriptome assembly and refinement (Trinity + PASA)
4. Open reading frame (ORF) prediction (TransDecoder)
5. Gene structure prediction using:
   - AUGUSTUS
   - GlimmerHMM
   - GeneMark-ET / GeneMark-EP+
6. Integration of multi-source evidence using EvidenceModeler (EVM)
7. Final gene model refinement, UTR annotation, and alternative splicing correction with PASA

A schematic representation of the HACGA1 workflow is shown below:

![HACGA1 workflow](figures/hacga1_pipeline.png)

---

## 📦 Deployment & Requirements

HACGA1 supports **Docker-based deployment only**.

### System Requirements
- Docker Engine
- Docker Compose

### Mandatory Requirement
- A valid **GeneMark license key**

> ⚠️ GeneMark is required for gene prediction. The pipeline will not run without a valid license.

---

## 🔑 GeneMark License Setup

1. Apply for a GeneMark license via the official [GeneMark website](https://genemark.bme.gatech.edu/GeneMark/license_download.cgi) or GeneMark@home.
2. Upon approval, save the license key to the following location:

```markdown
~/.gm_key
```

## 🚀 Installation

Clone the repository:

```shell
git clone https://github.com/liqianQ1/HACGA1.git
cd annotation_gene
```

Build the Docker image:

```shell
docker-compose build
```

------

## ▶️ Running HACGA1

Container runtime parameters and data volume mounts are defined in `docker-compose.yml`.
 Modify the mounted data directory to match your local data layout.

### Option 1: Run inside the container

```shell
docker exec -it hacga1 /bin/bash

python /pipeline/annotation_gene/ann_gene.py \
  --Outputdir /data \
  --Genome /data/tab/genome.tab \
  --Homolog /data/tab/homolog.tab \
  --RNAseq /data/tab/RNAseq.tab \
  --EST /data/tab/EST.tab \
  --max_parallel 100
```

### Option 2: Run directly from the host

```shell
docker exec hacga1 python /pipeline/annotation_gene/ann_gene.py \
  --Outputdir /data \
  --Genome /data/tab/genome.tab \
  --Homolog /data/tab/homolog.tab \
  --RNAseq /data/tab/RNAseq.tab \
  --EST /data/tab/EST.tab \
  --max_parallel 100
```

------

## 📥 Input & Output Data

HACGA1 requires the following categories of input:

- Genome sequence and repeat annotation
- RNA-seq reads
- Trinity-assembled transcript sequences
- Homologous protein sequences
- Configuration and evidence weight files

The pipeline generates structured outputs for each annotation stage, including:

- RNA-seq-based transcript annotations
- De novo gene predictions
- Homology-supported gene models
- Consensus gene structures integrated by EVM
- PASA-refined gene models with UTRs and alternative splicing


Detailed descriptions of input formats are provided in:

```
docs/data_format.md
```

------


## ⚠️ Important Notes

### Persistence of Trained AUGUSTUS Models

Trained AUGUSTUS models are stored inside the container and **will be removed when the container is deleted**.

- Default internal path:

```
/opt/conda/envs/ann/config/
```

- Recommended persistent storage location:

```
/data/augustus_models/
```

------

## 📬 Support & Contact

For questions, bug reports, or feature requests, please open an issue on GitHub.

------


**HACGA1** — *A scalable and reproducible solution for high-quality genome annotation.*



