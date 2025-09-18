#!/usr/bin/env python3
# targeted_brca.py — map → assemble → panel counts → plot 

import os, subprocess as sp, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = os.environ.get("BASE", "/content/drive/MyDrive/brca-targeted-analysis/brca_rna_targeted")
KEEP = os.environ.get("KEEP", "ERR13137440")

PROC = f"{BASE}/data/processed"
ALN  = f"{BASE}/align"
ISO  = f"{BASE}/isoforms"
RES  = f"{BASE}/results"
REFD = f"{BASE}/refs"

FA    = f"{REFD}/GRCh38.primary_assembly.genome.fa"
GTF   = f"{REFD}/gencode.v46.annotation.gtf"
PANEL = f"{REFD}/panel28_genes.bed"
FQXZ  = f"{PROC}/mini_{KEEP}.fastq.xz"

for d in (ALN, f"{ISO}/{KEEP}", RES):
    os.makedirs(d, exist_ok=True)

assert os.path.exists(FA),    f"missing {FA}"
assert os.path.exists(GTF),   f"missing {GTF}"
assert os.path.exists(PANEL), f"missing {PANEL}"
assert os.path.exists(FQXZ),  f"missing {FQXZ}"

def sh(cmd: str):
    print("+", cmd, flush=True)
    sp.run(["bash","-lc", cmd], check=True)

# 1) map .xz → BAM
sh(f"xzcat '{FQXZ}' | minimap2 -t 4 -ax splice -uf -k14 '{FA}' - | "
   f"samtools sort -m 1G -@2 -o '{ALN}/{KEEP}.bam' -")
sh(f"samtools index '{ALN}/{KEEP}.bam'")

# 2) assemble and make GTF.tmp
sh(f"stringtie -L -G '{GTF}' -o '{ISO}/{KEEP}/{KEEP}.gtf' '{ALN}/{KEEP}.bam'")
sh(f"stringtie -e -B -G '{GTF}' -o '{ISO}/{KEEP}/{KEEP}.guided.gtf' '{ALN}/{KEEP}.bam'")
sh(f"gffread -E '{ISO}/{KEEP}/{KEEP}.gtf' -T -o '{ISO}/{KEEP}/{KEEP}.gtf.tmp'")

# 3) mawk-safe transcripts.bed
gtf_tmp = f"{ISO}/{KEEP}/{KEEP}.gtf.tmp"
bed_tx  = f"{ISO}/{KEEP}/{KEEP}.transcripts.bed"
sh(
    r"awk 'BEGIN{FS=OFS=""\t""} $3==""transcript""{"
    r"tid=""NA""; if (match($9,/transcript_id ""[^""]+""/)) tid=substr($9,RSTART+15,RLENGTH-16);"
    r"print $1,$4-1,$5,tid,""."",$7"
    r"}' " + f"'{gtf_tmp}' > '{bed_tx}'"
)

# 4) counts and panel intersect
n_all = int(sp.check_output(["bash","-lc", f"grep -c $'\\ttranscript\\t' '{gtf_tmp}' || echo 0"]).decode().strip() or 0)
with open(f"{RES}/isoform_counts.tsv","w") as f:
    f.write(f"{KEEP}\t{n_all}\t0\n")

sh(f"bedtools intersect -s -wa -wb -a '{bed_tx}' -b '{PANEL}' > '{ISO}/{KEEP}/{KEEP}.transcripts.panel.bed' || true")
n_panel = int(sp.check_output(["bash","-lc", f"wc -l < '{ISO}/{KEEP}/{KEEP}.transcripts.panel.bed' 2>/dev/null || echo 0"]).decode().strip() or 0)
with open(f"{RES}/isoform_counts.tsv","w") as f:
    f.write(f"{KEEP}\t{n_all}\t{n_panel}\n")

# tidy and aggregate
rows = []
with open(f"{ISO}/{KEEP}/{KEEP}.transcripts.panel.bed") as fin:
    for ln in fin:
        if not ln.strip(): continue
        f = ln.rstrip("\n").split("\t")
        if len(f) < 12: continue
        rows.append((KEEP, f[9], 1))
tidy = pd.DataFrame(rows, columns=["sample","gene","count"])
tidy.to_csv(f"{RES}/panel_gene_hits.tsv", sep="\t", header=False, index=False)

agg = (tidy.groupby("gene", as_index=False)["count"]
            .sum().sort_values("count", ascending=False))
agg.to_csv(f"{RES}/panel_gene_hits_agg.tsv", sep="\t", index=False)

# 5) plot
plt.figure(figsize=(8,5))
agg.head(15).set_index("gene")["count"][::-1].plot(kind="barh")
plt.xlabel("Transcript overlaps (panel)")
plt.title("Top HBOC panel genes (subsampled)")
plt.tight_layout()
plt.savefig(f"{RES}/top_panel_genes.png", dpi=160, bbox_inches="tight")

print("[DONE]")
print(open(f"{RES}/isoform_counts.tsv").read().strip())
