# BRCA-Targeted RNA Analysis: HBOC 28-Gene Panel (ONT Long-Read) mini analysis

## 1. Data sources
This workflow reproduces part of the SOSTAR pipeline using a publicly available capture-based targeted RNA-seq dataset focused on 28 hereditary breast/ovarian cancer (HBOC) genes. 
The sequencing was performed on ONT MinION runs, and we use a subsampled FASTQ (<1 GB) to ensure reproducibility in Colab.

**Paper:** Fine mapping of RNA isoform diversity using an innovative targeted long-read RNA sequencing protocol with novel dedicated bioinformatics pipeline. 
BMC Genomics, 2024. **DOI:** [10.1186/s12864-024-10741-0](https://doi.org/10.1186/s12864-024-10741-0)

**Study accession:** ENA PRJEB75906 (28 HBOC genes, lymphoblastoid cell lines);
**Download:** 
Example run ERR13137440 FASTQ from ENA : 
	•	HTTPS: https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR131/040/ERR13137440/ERR13137440.fastq.gz
	•	FTP: ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR131/040/ERR13137440/ERR13137440.fastq.gz

### Layout
```bash
brca-targeted-analysis/
├── README.md
└── brca-targeted-analysis/
    ├── data/
    │   ├── raw/                # original FASTQ(s) from ENA/SRA
    │   └── processed/          # subsampled FASTQ(s) (<1 GB each)
    ├── refs/                   # GRCh38 + GENCODE + 28-gene BED
    ├── align/                  # BAM/CRAM alignments
    ├── isoforms/               # StringTie assemblies + transcript BEDs
    ├── results/                # tables + plots
    │   ├── isoform_counts.tsv
    │   ├── panel_gene_hits.tsv
    │   └── top_panel_genes.png
    ├── metadata.yaml
    ├── questions.yaml
    ├── answers.yaml
    └── workflow/
        ├──targeted_brca.py     # python file
        └── outputs/            # extra figs/zips 
```

## 2. Installation & Execution
### Option A: Run in Google Colab (recommended)
	1.	Open a Colab session.
	2.	Mount Google Drive:
 ```python
from google.colab import drive
drive.mount('/content/drive')
```
	3.	Install required packages:
 ```bash
!apt-get -qq update
!apt-get -qq install -y minimap2 samtools bedtools seqtk gffread stringtie
```
	4.	Place raw FASTQ in data/raw/ and subsample:
 ```bash
!seqtk sample -s100 data/raw/ERR13137440.fastq.gz 300000 | gzip -c > data/processed/mini_ERR13137440.fastq.gz
```
	5.	Run the workflow:
 ```bash
!python brca-targeted-analysis/workflow/targeted_brca.py
```
#### Option B: Run locally / HPC
	1.	Clone repo:
 ```bash
git clone https://github.com/<your-username>/brca-targeted-analysis.git
cd brca-targeted-analysis/brca-targeted-analysis/workflow
```
	2.	Install dependencies (apt or conda):
	•	minimap2
	•	samtools
	•	bedtools
	•	seqtk
	•	gffread
	•	stringtie
	•	python ≥3.10 with pandas, matplotlib, pyyaml

 	3.	Provide raw FASTQ under data/raw/, subsample with seqtk.
	4.	Run:
 ```bash
python targeted_brca.py
```
 ## 3. Primary Outputs
 	•	results/isoform_counts.tsv → all transcripts + panel overlaps.
	•	results/panel_gene_hits.tsv → per-gene transcript overlap counts.
	•	results/top_panel_genes.png → barplot of top HBOC genes (e.g. ATM, BRCA1, RAD50).

All outputs are under results/, all intermediates under align/ and isoforms/.

## Notes:
	•	Input FASTQs are subsampled to 300k reads to keep files <1 GB.
	•	Pipeline uses only standard tools (minimap2, samtools, stringtie, gffread, bedtools).
	•	Entire workflow is captured in a single Python script (targeted_brca.py)
	•	Questions and answers are generated automatically (questions.yaml, answers.yaml).

