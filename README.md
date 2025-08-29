# BRCA-Targeted RNA Analysis: HBOC 28-Gene Panel (ONT Long-Read) mini analysis

## 1. Data sources
This workflow reproduces part of the SOSTAR pipeline using a publicly available capture-based targeted RNA-seq dataset focused on 28 hereditary breast/ovarian cancer (HBOC) genes. 
The sequencing was performed on ONT MinION runs, and we use a subsampled FASTQ (<1 GB) to ensure reproducibility in Colab.

**Paper:** Fine mapping of RNA isoform diversity using an innovative targeted long-read RNA sequencing protocol with novel dedicated bioinformatics pipeline. 
BMC Genomics, 2024. **DOI:** [10.1186/s12864-024-10741-0](https://doi.org/10.1186/s12864-024-10741-0)

**Study accession:** ENA PRJEB75906 (28 HBOC genes, lymphoblastoid cell lines);
**Download:** Example run ERR13137440 FASTQ from ENA : ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR131/040/ERR13137440/ERR13137440.fastq.gz 

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
        ├── targeted_brca.ipynb      # Colab notebook (end-to-end)
        └── outputs/            # extra figs/zips
```

## 2. Installation & Execution
### Step 1. Open workflow/targeted_brca.ipynb in Google Colab.
### Step 2. Mount Google Drive:
```python
from google.colab import drive
drive.mount('/content/drive')
```
### Step 3. Install required packages (done inside notebook):
```bash
apt-get -qq update
apt-get -qq install -y minimap2 samtools bedtools seqtk gffread stringtie
```
### Step 4. Place raw FASTQ in data/raw/ and subsample:
```bash
seqtk sample -s100 data/raw/ERR13137440.fastq.gz 300000 | gzip -c > data/processed/mini_ERR13137440.fastq.gz
```
### Step 5. Run the notebook cells to:
	•	Build refs (GRCh38 + GENCODE v46 + 28-gene BED).
	•	Map reads → BAM (minimap2 + samtools).
	•	Assemble isoforms (StringTie).
	•	Convert to BED (gffread + Python).
	•	Intersect with panel BED (bedtools).
	•	Summarize and plot results.

 ## 3. Primary Outputs
 	•	results/isoform_counts.tsv → all transcripts + panel overlaps.
	•	results/panel_gene_hits.tsv → per-gene transcript overlap counts.
	•	results/top_panel_genes.png → barplot of top HBOC genes (e.g. ATM, BRCA1, RAD50).

All outputs are under results/, all intermediates under align/ and isoforms/.

## Notes:
	•	Input FASTQs are subsampled to 300k reads to keep files <1 GB.
	•	Pipeline uses only standard tools: minimap2, samtools, stringtie, gffread, bedtools.
	•	Reproducibility: all steps contained in Colab notebook workflow.ipynb.
	•	You can scale to 1M reads or add more runs if resources allow.

