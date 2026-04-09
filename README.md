# HLA-I Multiomics Analysis in HCC

![header](doc/header.png)

## 📄 Associated Manuscript

This repository contains the analysis code for the study:

**"HLA Class I Expression and Ratios as Diagnostic and Predictive Biomarkers in Hepatocellular Carcinoma: A Multi-Omics Study"**

🟡 Currently under publication in *Heliyon*.

## 📁 Project Structure

```
HLA-I-Multiomics-HCC/
├── README.md
├── LICENSE
├── doc/
│   └── header.png
├── data/
│   ├── bulk_RNA_seq/
│   ├── differential_correlation/
│   ├── genotype_allele_specific_expression/
│   ├── model/
│   ├── mutation/
│   ├── neoantigen/
│   ├── scrna/
│   └── spatial_transcriptomics/
├── results/
└── scripts/
```

---

## ⚠️ Data Availability

Due to file size and licensing restrictions, datasets are NOT included.

Users must download them manually:

### Bulk expression:
https://xenabrowser.net/datapages/?dataset=TcgaTargetGtex_rsem_gene_tpm&host=https://toil.xenahubs.net

### scRNA-seq & spatial:
http://lifeome.net:809/#/download

### Processed data

Some of the datasets used in this study, whose original references are cited in the manuscript, were further cleaned, matched with metadata, and subsetted for HLA-focused analyses. These processed files are available at the following link:

https://drive.google.com/file/d/1xofyMsXwu7iNUZW-dApgG6b46AvYE42Y/view?usp=sharing

After downloading and unzipping the archive, the extracted files should be placed in their corresponding folders within the `data/` directory of this repository.

---

## ▶️ How to Run

```r
setwd("path/to/HLA-I-Multiomics-HCC")
source("scripts/01_allele_frequency.R")
```

---

## 🧬 Pipeline Overview

1. Allele frequency
2. Allele expression
3. Neoantigen affinity
4. Clinical analysis
5. Correlation analysis
6. Machine learning
7. Mutation analysis
8. Single-cell & spatial

---

## 🧪 Neoantigen Analysis (netMHCpan)

Neoantigen binding affinity prediction is performed using **netMHCpan**.  
The pipeline automatically generates peptide files and corresponding shell scripts for each HLA locus (A, B, C).

### 🔹 Step 1: Prepare input files

After running:

```r
source("scripts/03_neoantigen_affinity.R")
```
the following will be generated:

Peptide files (*.pep.text)
Shell scripts:
Code-Linux-HLA-A.sh
Code-Linux-HLA-B.sh
Code-Linux-HLA-C.sh

These are typically located in:

results/neoantigen/linux_scripts/
results/neoantigen/peptides/
🔹 Step 2: Run netMHCpan (Linux required)

Transfer the peptide files and shell scripts to a Linux environment where netMHCpan is installed.

Run:

bash Code-Linux-HLA-A.sh
bash Code-Linux-HLA-B.sh
bash Code-Linux-HLA-C.sh

Each script will:

Loop over samples/alleles
Run netMHCpan for each peptide set
Generate raw binding affinity outputs
🔹 Step 3: Collect output files

After execution, you will obtain multiple output files (typically .txt or .xls depending on netMHCpan version).

These files should be organized into the following directories:

results/neoantigen/A/
results/neoantigen/B/
results/neoantigen/C/

Each folder must contain outputs corresponding to that locus.

🔹 Step 4: Downstream analysis in R

Once outputs are placed correctly, re-run:

source("scripts/03_neoantigen_affinity.R")

This step will:

Parse netMHCpan output files
Extract binding affinity values
Merge with HLA genotype and peptide data
Perform statistical analysis:
ANOVA
Post-hoc tests

Summary statistics
⚠️ Important Notes
netMHCpan is not included in this repository and must be installed separately.
Ensure that:
File names are not modified
Outputs are placed in the correct directories
The pipeline assumes consistent naming between peptide files and output files.


## 📚 Citation

If you use this repository, please cite:

Aliyari S. et al. *Heliyon* (in press)

---

## 🧑‍💻 Author

Shahram Aliyari  
Bioinformatics Researcher
