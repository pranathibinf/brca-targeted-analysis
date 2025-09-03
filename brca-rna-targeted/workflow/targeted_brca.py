#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# targeted_brca.py 

import os, sys, re, shutil, gzip, urllib.request, subprocess as sp
from pathlib import Path

def sh(cmd, **kw):
    print("+", cmd if isinstance(cmd, str) else " ".join(cmd), flush=True)
    return sp.run(cmd if isinstance(cmd, list) else ["bash","-lc",cmd], check=True, **kw)

# --- Mount Drive (if needed and running on colab) ---
def maybe_mount_colab_drive():
    try:
        import google.colab  # type: ignore
        from google.colab import drive  # type: ignore
        drive.mount('/content/drive')
        return True
    except Exception:
        return False

IN_COLAB = maybe_mount_colab_drive()

# --- Setup BASE/KEEP and directories ---
BASE = "/content/drive/MyDrive/brca-targeted-analysis/brca_rna_targeted" if IN_COLAB \
       else os.environ.get("BASE", os.path.abspath("./brca-targeted-analysis/brca_rna_targeted"))
KEEP = "ERR13137440"  # single sample

for d in ["data/raw","data/processed","refs","align","isoforms","results","workflow/outputs"]:
    os.makedirs(os.path.join(BASE, d), exist_ok=True)

os.environ["BASE"] = BASE
os.environ["KEEP"] = KEEP

print(BASE)
for root, dirs, files in os.walk(BASE):
    level = root.replace(BASE, "").count(os.sep)
    indent = " " * (2 * level)
    print(f"{indent}{os.path.basename(root)}/")
    subindent = " " * (2 * (level + 1))
    for f in files:
        print(f"{subindent}{f}")

# --- (INSTALL TOOLS) disabled by default ---
# To run, set RUN_INSTALL=True
RUN_INSTALL = False
if RUN_INSTALL:
    sh("apt-get -qq update")
    sh("apt-get -qq install -y minimap2 samtools bedtools seqtk gffread stringtie")

# --- (SUBSAMPLE)  ---
RUN_SUBSAMPLE = False  # set True if you want to run subsampling here
if RUN_SUBSAMPLE:
    cmd = f'''
BASE="{BASE}"; KEEP="{KEEP}"
RAW="$BASE/data/raw"; PROC="$BASE/data/processed"; NREADS=300000
if [ -s "$RAW/${{KEEP}}.fastq.gz" ]; then
  seqtk sample -s100 "$RAW/${{KEEP}}.fastq.gz" "$NREADS" | gzip -c > "$PROC/mini_${{KEEP}}.fastq.gz"
fi
ls -lh "$PROC" || true
'''
    sh(cmd)

# --- Build refs (GRCh38, GENCODE v46) and 28-gene BED ---
REF = f"{BASE}/refs"
os.makedirs(REF, exist_ok=True)
fa_gz="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/GRCh38.primary_assembly.genome.fa.gz"
gtf_gz="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.annotation.gtf.gz"
fa=f"{REF}/GRCh38.primary_assembly.genome.fa"
gtf=f"{REF}/gencode.v46.annotation.gtf"

def fetch_gz(url, out_plain):
    if not os.path.exists(out_plain):
        gz = out_plain + ".gz"
        urllib.request.urlretrieve(url, gz)
        with gzip.open(gz,'rb') as fin, open(out_plain,'wb') as fout: shutil.copyfileobj(fin,fout)

fetch_gz(fa_gz, fa)
fetch_gz(gtf_gz, gtf)

panel_txt=f"{REF}/panel28_genes.txt"
with open(panel_txt,"w") as f:
    f.write("\n".join([
      "BRCA1","BRCA2","ATM","ATR","BARD1","BRIP1","CDH1","CHEK2","MRE11","NBN","PALB2","PTEN",
      "RAD50","RAD51","RAD51C","RAD51D","RECQL","TP53","XRCC2","XRCC3","FANCA","FANCC","FANCM",
      "MSH2","MSH6","MLH1","PMS2","EPCAM"
    ])+"\n")

raw_bed=f"{REF}/panel28_genes.raw.bed"
bed_bed=f"{REF}/panel28_genes.bed"
genes=set(x.strip() for x in open(panel_txt) if x.strip())
n=0
with open(gtf) as fin, open(raw_bed,"w") as fout:
    for line in fin:
        if not line or line[0]=="#": continue
        parts=line.rstrip("\n").split("\t")
        if len(parts)<9 or parts[2]!="gene": continue
        m=re.search(r'gene_name "([^"]+)"', parts[8])
        if not m: continue
        name=m.group(1)
        if name in genes:
            chrom=parts[0]; start=int(parts[3])-1; end=int(parts[4]); strand=parts[6]
            fout.write(f"{chrom}\t{start}\t{end}\t{name}\t.\t{strand}\n")
            n+=1
shutil.copy(raw_bed, bed_bed)
print("panel28 genes BED rows:", n)

# --- Map (minimap2) ---
RUN_MAP = True
if RUN_MAP:
    cmd = f'''
BASE="{BASE}"; KEEP="{KEEP}"
PROC="$BASE/data/processed"; ALN="$BASE/align"; REF="$BASE/refs"
FA="$REF/GRCh38.primary_assembly.genome.fa"

if [ -s "$ALN/${{KEEP}}.bam" ]; then
  echo "[INFO] BAM exists: $ALN/${{KEEP}}.bam"
else
  fq="$PROC/mini_${{KEEP}}.fastq.gz"
  if [ -s "$fq" ]; then
    minimap2 -t 4 -ax splice -uf -k14 "$FA" "$fq" | samtools sort -m 1G -@2 -o "$ALN/${{KEEP}}.bam" -
    samtools index "$ALN/${{KEEP}}.bam"
  else
    echo "[WARN] Missing $fq — provide subsampled FASTQ or skip to annotation if BAM exists."
  fi
fi
ls -lh "$ALN" || true
'''
    sh(cmd)

# --- Annotate (StringTie) and make transcript BED  ---
RUN_ANNOTATE = True
if RUN_ANNOTATE:
    cmd = f'''
BASE="{BASE}"; KEEP="{KEEP}"
ALN="$BASE/align"; ISO="$BASE/isoforms"; RES="$BASE/results"; REF="$BASE/refs"
GTF="$REF/gencode.v46.annotation.gtf"

bam="$ALN/${{KEEP}}.bam"; mkdir -p "$ISO/$KEEP" "$RES"
[ -s "$bam" ] || {{ echo "[ERR] BAM missing: $bam"; exit 1; }}

stringtie -L -G "$GTF" -o "$ISO/$KEEP/${{KEEP}}.gtf" "$bam"
stringtie -e -B -G "$GTF" -o "$ISO/$KEEP/${{KEEP}}.guided.gtf" "$bam"
gffread -E "$ISO/$KEEP/${{KEEP}}.gtf" -T -o "$ISO/$KEEP/${{KEEP}}.gtf.tmp"
'''
    sh(cmd)

    # inline python (as in notebook) to write transcripts.bed
    gtf_tmp=f"{BASE}/isoforms/{KEEP}/{KEEP}.gtf.tmp"
    bed_out=f"{BASE}/isoforms/{KEEP}/{KEEP}.transcripts.bed"
    n=0
    with open(gtf_tmp) as fin, open(bed_out,"w") as fout:
        for ln in fin:
            if not ln or ln[0]=="#": continue
            f=ln.rstrip("\n").split("\t")
            if len(f)<9 or f[2]!="transcript": continue
            m=re.search(r'transcript_id "([^"]+)"', f[8])
            tid=m.group(1) if m else "NA"
            chrom=f[0]; start=str(int(f[3])-1); end=f[4]; strand=f[6]
            fout.write("\t".join([chrom,start,end,tid,".",strand])+"\n")
            n+=1
    print("transcript BED rows:", n)

    # quick global counts 
    iso_counts=f"{BASE}/results/isoform_counts.tsv"
    open(iso_counts,"w").close()
    # count transcript lines
    n_all = 0
    with open(gtf_tmp) as f:
        for ln in f:
            if "\ttranscript\t" in ln:
                n_all += 1
    with open(iso_counts,"w") as f:
        f.write(f"{KEEP}\t{n_all}\t0\n")
    print(open(iso_counts).read())

# --- Intersect with 28-gene panel BED & per-gene counts  ---
RUN_INTERSECT = True
if RUN_INTERSECT:
    cmd = f'''
BASE="{BASE}"; KEEP="{KEEP}"
ISO="$BASE/isoforms"; RES="$BASE/results"; REF="$BASE/refs"
PANELBED="$REF/panel28_genes.bed"

bedtools intersect -s -wa -wb -a "$ISO/$KEEP/${{KEEP}}.transcripts.bed" -b "$PANELBED" \
  > "$ISO/$KEEP/${{KEEP}}.transcripts.panel.bed"

n_panel=$(wc -l < "$ISO/$KEEP/${{KEEP}}.transcripts.panel.bed" || echo 0)
awk -v OFS="\\t" -v S="$KEEP" -v P="$n_panel" 'NR==1{{print S,$2,P}}' "$RES/isoform_counts.tsv" > "$RES/.tmp" || echo -e "$KEEP\\t0\\t$P" > "$RES/.tmp"
mv "$RES/.tmp" "$RES/isoform_counts.tsv"

: > "$RES/panel_gene_hits.tsv"
bedtools intersect -s -wa -wb -a "$ISO/$KEEP/${{KEEP}}.transcripts.bed" -b "$PANELBED" \
 | awk -v S="$KEEP" 'BEGIN{{FS=OFS="\\t"}}{{print S,$10,$4}}' \
 | sort | uniq -c \
 | awk 'BEGIN{{OFS="\\t"}}{{print $2,$3,$1}}' >> "$RES/panel_gene_hits.tsv"

echo "[RESULTS]"; cat "$RES/isoform_counts.tsv"
echo "[GENE-HIT HEAD]"; head -n 15 "$RES/panel_gene_hits.tsv" || true
'''
    sh(cmd)

# --- Plot top genes and save PNG ---
import pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

tsv=os.path.join(BASE,"results","panel_gene_hits.tsv")
df=pd.read_csv(tsv, sep="\t", header=None, names=["sample","gene","count"], dtype=str)
df["count"]=pd.to_numeric(df["count"], errors="coerce")
df=df.dropna(subset=["count"])
summary=df.groupby("gene", as_index=False)["count"].sum().sort_values("count", ascending=False)

print("Top 10 genes:\n", summary.head(10).to_string(index=False))

top15=summary.head(15).set_index("gene")
plt.figure(figsize=(8,5))
top15["count"].plot(kind="barh")
plt.xlabel("Transcript overlaps (panel)")
plt.title("Top 15 HBOC Panel Genes (subsampled)")
plt.tight_layout()
out=os.path.join(BASE,"results","top_panel_genes.png")
plt.savefig(out, dpi=160, bbox_inches="tight")
print("saved:", out)

# --- (OPTIONAL) BAM→CRAM ---
RUN_CRAM = False
if RUN_CRAM:
    cmd = f'''
BASE="{BASE}"; KEEP="{KEEP}"
ALN="$BASE/align"; REF="$BASE/refs/GRCh38.primary_assembly.genome.fa"
bam="$ALN/${{KEEP}}.bam"
if [ -s "$bam" ]; then
  cram="${{bam%.bam}}.cram"
  samtools view -T "$REF" -C -o "$cram" "$bam"
  samtools index "$cram"
  rm -f "$bam" "$bam.bai"
fi
ls -lh "$ALN" || true
'''
    sh(cmd)

# --- write questions.yaml ---
import textwrap
qs_path = os.path.join(BASE, "questions.yaml")
questions_yaml = textwrap.dedent("""\
task: >
  You have been given a subsampled long-read RNA-seq dataset (ERR13137440).
  Using the processed FASTQ/BAM and panel-intersect outputs in results/, answer the following:

questions:
- id: q1_total_isoforms
  stage: isoform_assembly
  text: How many total transcripts were assembled from ERR13137440?
  answer_type: integer_exact

- id: q2_panel_isoforms
  stage: panel_intersection
  text: How many transcripts overlapped the 28-gene HBOC panel?
  answer_type: integer_exact

- id: q3_top_gene
  stage: gene_summary
  text: Which panel gene had the highest transcript overlap count?
  answer_type: string_exact

- id: q4_brca1_vs_brca2
  stage: gene_comparison
  text: What is the transcript count ratio BRCA1/BRCA2?
  answer_type: numeric_float
  tolerance: 0.01

- id: q5_top3_genes
  stage: visualization
  text: According to panel_gene_hits.tsv, what are the top 3 genes by count?
  answer_type: list_string
""")
with open(qs_path, "w") as f:
    f.write(questions_yaml)
print(f"[OK] wrote {qs_path}")
print(open(qs_path).read())

# --- compute answers.yaml  ---
import math, yaml
iso_tsv = os.path.join(BASE, "results", "isoform_counts.tsv")
hits_tsv = os.path.join(BASE, "results", "panel_gene_hits.tsv")
ans_path = os.path.join(BASE, "answers.yaml")

iso = pd.read_csv(iso_tsv, sep="\t", header=None, names=["sample","total","panel"])
q1_total_isoforms = int(iso.loc[0, "total"])
q2_panel_isoforms = int(iso.loc[0, "panel"])

hits = pd.read_csv(hits_tsv, sep="\t", header=None, names=["sample","gene","count"])
hits["count"] = pd.to_numeric(hits["count"], errors="coerce").fillna(0).astype(int)
gene_sum = hits.groupby("gene", as_index=True)["count"].sum().sort_values(ascending=False)

q3_top_gene = str(gene_sum.index[0]) if not gene_sum.empty else ""

brca1 = int(gene_sum.get("BRCA1", 0))
brca2 = int(gene_sum.get("BRCA2", 0))
q4_brca1_vs_brca2 = float(brca1)/float(brca2) if brca2 != 0 else float("inf")

top3 = list(gene_sum.head(3).index)
answers = {
    "q1_total_isoforms": int(q1_total_isoforms),
    "q2_panel_isoforms": int(q2_panel_isoforms),
    "q3_top_gene": q3_top_gene,
    "q4_brca1_vs_brca2": None if math.isinf(q4_brca1_vs_brca2) else round(q4_brca1_vs_brca2, 2),
    "q5_top3_genes": top3,
}
with open(ans_path, "w") as f:
    yaml.safe_dump(answers, f, sort_keys=False)
print(f"[OK] wrote {ans_path}")
print(open(ans_path).read())
