# HLA-I Multiomics Analysis in HCC

![header](doc/header.png)

## 📄 Associated Manuscript

This repository contains the analysis code for the study:

**"HLA Class I Expression and Ratios as Diagnostic and Predictive Biomarkers in Hepatocellular Carcinoma: A Multi-Omics Study"**

🟡 Currently under publication in *Heliyon*.

---

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

Datasets are not included due to size and licensing restrictions.

### Download links:

- Bulk expression (TCGA + GTEx):
https://xenabrowser.net/datapages/?dataset=TcgaTargetGtex_rsem_gene_tpm&host=https://toil.xenahubs.net

- scRNA-seq and spatial transcriptomics:
http://lifeome.net:809/#/download

---

### 🔬 Processed Data

Some datasets used in this study were:

- cleaned  
- matched with metadata  
- filtered for HLA-focused analyses  

Download processed data here:
https://drive.google.com/file/d/1xofyMsXwu7iNUZW-dApgG6b46AvYE42Y/view?usp=sharing

After extraction, place files in the appropriate folders inside `data/`.

---

## ▶️ How to Run

Always run scripts from the project root:

```r
setwd("path/to/HLA-I-Multiomics-HCC")
source("scripts/01_allele_frequency.R")
```

---

## 🧬 Pipeline Overview

| Step | Script | Description |
|------|--------|-------------|
| 01 | allele_frequency | Calculate HLA allele frequencies |
| 02 | allele_expression | Compare allele-specific expression |
| 03 | neoantigen_affinity | Predict peptide binding affinity |
| 04 | HCCDB_expression | Expression comparison (HCC vs normal) |
| 05 | clinical_analysis | Clinical and survival analysis |
| 06 | correlation | GTEx vs TCGA correlation |
| 07 | biomarker | Diagnostic model (logistic regression) |
| 08 | mutation | Mutation landscape (TCGA) |
| 09 | scRNA | Single-cell analysis |
| 10 | spatial | Spatial transcriptomics |

---

## 🧪 Neoantigen Analysis (netMHCpan)

### Step 1: Generate inputs

```r
source("scripts/03_neoantigen_affinity.R")
```

Generates:
- peptide files (`*.pep.text`)
- shell scripts

### Step 2: Run in Linux

```bash
bash Code-Linux-HLA-A.sh
bash Code-Linux-HLA-B.sh
bash Code-Linux-HLA-C.sh
```

### Step 3: Move outputs

```
results/neoantigen/A/
results/neoantigen/B/
results/neoantigen/C/
```

### Step 4: Parse results

```r
source("scripts/03_neoantigen_affinity.R")
```

---

## 📚 Citation

Aliyari S. et al. *Heliyon* (in press)

---

## 🧑‍💻 Author

Shahram Aliyari  
Bioinformatics Researcher
