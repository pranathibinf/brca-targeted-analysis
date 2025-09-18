# BRCA-Targeted RNA Analysis: HBOC 28-Gene Panel (ONT Long-Read) mini analysis

## 1. Data sources
This workflow reproduces part of the SOSTAR pipeline using a publicly available capture-based targeted RNA-seq dataset focused on 28 hereditary breast/ovarian cancer (HBOC) genes.  
Sequencing was performed on ONT MinION runs. We use a **subsampled FASTQ.xz (~50k–300k reads, <200 MB)** for reproducibility.

**Paper:** Fine mapping of RNA isoform diversity using an innovative targeted long-read RNA sequencing protocol with novel dedicated bioinformatics pipeline.  
*BMC Genomics*, 2024. DOI: [10.1186/s12864-024-10741-0](https://doi.org/10.1186/s12864-024-10741-0)

**Study accession:** ENA PRJEB75906  
**Example run:** ERR13137440 FASTQ

- HTTPS: https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR131/040/ERR13137440/ERR13137440.fastq.gz  
- FTP: ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR131/040/ERR13137440/ERR13137440.fastq.gz  

## 2. Layout

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
---

## 3. Installation & Execution

### Option A: Google Colab (testing only)
1. Open a Colab session.  
2. Mount Google Drive:
```python
from google.colab import drive
drive.mount('/content/drive')
3. Place raw FASTQ under data/raw/.
Use the provided script fetch_refs_and_data.sh or the notebook snippet to subsample and convert to .xz:
 ```bash
!bash brca-targeted-analysis/workflow/fetch_refs_and_data.sh
```
4.	Run the pipeline:
```bash
!python brca-targeted-analysis/workflow/targeted_brca.py
```

#### Option B: Run locally / HPC
	1.	Clone repo:
 ```bash
git clone https://github.com/pranathibinf/brca-targeted-analysis.git
cd brca-targeted-analysis/brca-targeted-analysis/workflow
```
	2.	Install dependencies:
	•	minimap2
	•	samtools
	•	bedtools
	•	stringtie
	•	gffread
	•	python ≥3.10 with pandas, matplotlib, pyyaml

 	3.	Provide raw FASTQ under data/raw/, and subsample with:
  ```bash
seqtk sample -s100 data/raw/ERR13137440.fastq.gz 50000 | gzip -c > data/processed/mini_ERR13137440.fastq.gz
python3 compress_to_xz.py  # converts .gz to .xz using lzma
```  
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
	•	Subsampled to ~50k–300k reads to reduce runtime and keep processed FASTQ <200 MB.
	•	Uses only standard tools.
	•	All steps are contained in a single Python script: workflow/targeted_brca.py.

