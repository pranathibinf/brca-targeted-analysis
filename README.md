# BRCA-Targeted RNA Analysis: HBOC 28-Gene Panel (ONT Long-Read) mini analysis

## 1. Data sources
This workflow reproduces part of the SOSTAR pipeline using a publicly available capture-based targeted RNA-seq dataset focused on 28 hereditary breast/ovarian cancer (HBOC) genes.  
Sequencing was performed on ONT MinION runs. We use a **subsampled FASTQ.gz (~15k reads)** for reproducibility.

**Paper:** Fine mapping of RNA isoform diversity using an innovative targeted long-read RNA sequencing protocol with novel dedicated bioinformatics pipeline.  
*BMC Genomics*, 2024. DOI: [10.1186/s12864-024-10741-0](https://doi.org/10.1186/s12864-024-10741-0)

**Study accession:** ENA PRJEB75906  
**Example run:** ERR13137440 FASTQ

- HTTPS: https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR131/040/ERR13137440/ERR13137440.fastq.gz  
- FTP: ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR131/040/ERR13137440/ERR13137440.fastq.gz  

For Taiga:
**Bundled references (panel-only, no downloads):**
- `refs/GRCh38.panel28.pad50k.fa`
- `refs/gencode.v46.panel28.gtf`
- `refs/panel28_genes.bed`

**Processed input provided:**
- `data/processed/mini_ERR13137440.fastq.gz`

## 2. Layout
```bash
brca-targeted-analysis/
├── README.md
└── brca-targeted-analysis/
    ├── data/
    │   ├── raw/                         # original FASTQs
    │   └── processed/
    │       └── mini_ERR13137440.fastq.gz (miniatured fastq)
    ├── refs/                            # bundled panel references
    │   ├── GRCh38.panel28.pad50k.fa
    │   ├── gencode.v46.panel28.gtf
    │   └── panel28_genes.bed
    ├── align/
    ├── isoforms/
    ├── results/
    │   ├── isoform_counts.tsv
    │   ├── panel_gene_hits.tsv
    │   └── top_panel_genes.png
    ├── metadata.yaml
    ├── questions.yaml
    ├── answers.yaml
    └── workflow/
        └── targeted_brca.py (workflow script)
```
## 3. Dependencies
- minimap2
- samtools
- stringtie
- gffread
- bedtools
- python ≥3.10 (pandas, matplotlib, pyyaml)
  
## 4. Installation & Execution
### Option A — Local / HPC (recommended for running)
```bash
# create env (conda based)
conda create -n brca-rna python=3.11 -y
conda activate brca-rna
conda config --add channels conda-forge
conda config --add channels bioconda
conda install -y minimap2 samtools bedtools stringtie gffread pandas matplotlib pyyaml
```
# run the workflow (set BASE to your local project root)
```bash
export BASE="/path/to/brca-targeted-analysis/brca-rna-targeted"
export KEEP="ERR13137440"
python "${BASE}/workflow/targeted_brca.py"
```
**Expected inputs present before running:**
	•	${BASE}/data/processed/mini_ERR13137440.fastq.xz
	•	${BASE}/refs/GRCh38.panel28.pad50k.fa
	•	${BASE}/refs/gencode.v46.panel28.gtf
	•	${BASE}/refs/panel28_genes.bed

**Outputs written to:**
	•	${BASE}/results/isoform_counts.tsv
	•	${BASE}/results/panel_gene_hits.tsv
	•	${BASE}/results/top_panel_genes.png
 
### Option B: Google colab (testing only)
```bash
from google.colab import drive
drive.mount('/content/drive')

BASE = "/content/drive/MyDrive/brca-targeted-analysis/brca_rna_targeted"
KEEP = "ERR13137440"

# run
!BASE="{BASE}" KEEP="{KEEP}" python "{BASE}/workflow/targeted_brca.py"
```
**Note:** The workflow uses the bundled panel references and a pre-subsampled .xz FASTQ. No downloads are performed by targeted_brca.py. If you want to regenerate the .xz from a full .fastq.gz, do that separately.

# 5. Primary Outputs
- `results/isoform_counts.tsv`
- `results/panel_gene_hits.tsv`
- `results/top_panel_genes.png`
  
## Notes:
	•	Subsampled to 15k reads to reduce runtime and keep processed FASTQ.
	•	Uses only standard tools.
	•	Workflow is captured in a single Python script: workflow/targeted_brca.py.

