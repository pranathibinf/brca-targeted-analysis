#!/usr/bin/env python3
import os, re, subprocess as sp, pandas as pd, matplotlib
from pathlib import Path
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def sh(cmd: str):
    import os, subprocess as sp
    env_path = os.environ.get("PATH", "")
    sp.run(["bash","-c", f'export PATH="{env_path}"; {cmd}'], check=True)

BASE = os.environ.get("BASE") or str(Path(__file__).resolve().parents[1])
KEEP = os.environ.get("KEEP", "ERR13137440")
ENV_BIN = os.environ.get("ENV_BIN", "")
if ENV_BIN:
    os.environ["PATH"] = ENV_BIN + os.pathsep + os.environ["PATH"]

REF  = f"{BASE}/refs"
PROC = f"{BASE}/data/processed"
ALN  = f"{BASE}/align"
ISOd = f"{BASE}/isoforms"
ISO  = f"{ISOd}/{KEEP}"
RES  = f"{BASE}/results"

for d in (ALN, ISO, RES): Path(d).mkdir(parents=True, exist_ok=True)

FA        = f"{REF}/GRCh38.panel28.pad50k.fa"
GTF       = f"{REF}/gencode.v46.panel28.gtf"
PANEL_G   = f"{REF}/panel28_genes.bed"
PANEL_P   = f"{REF}/panel28_genes.panelspace.bed"
FQGZ      = f"{PROC}/mini_{KEEP}.fastq.gz"
BAM       = f"{ALN}/{KEEP}.bam"
GTF_OUT   = f"{ISO}/{KEEP}.gtf"
GTF_TMP   = f"{ISO}/{KEEP}.gtf.tmp"
BED_TX    = f"{ISO}/{KEEP}.transcripts.bed"
BED_INT   = f"{ISO}/{KEEP}.transcripts.panel.bed"
PNG_OUT   = f"{RES}/top_panel_genes.png"
TSV_COUNTS= f"{RES}/isoform_counts.tsv"
TSV_HITS  = f"{RES}/panel_gene_hits.tsv"
TSV_AGG   = f"{RES}/panel_gene_hits_agg.tsv"

assert os.path.exists(FA),   f"missing {FA}"
assert os.path.exists(GTF),  f"missing {GTF}"
assert os.path.exists(PANEL_G), f"missing {PANEL_G}"
assert os.path.exists(FQGZ), f"missing {FQGZ}"

if not os.path.exists(PANEL_P):
    gene2info = {}
    with open(FA) as f:
        for ln in f:
            if not ln.startswith(">"): continue
            h=ln[1:].strip()
            m=re.match(r'^([^|]+)\|([^:]+):(\d+)-(\d+)\(([+-])\)', h)
            if not m: continue
            gene=m.group(1); gstart=int(m.group(3)); strand=m.group(5)
            contig=h
            gene2info[gene]=(contig, gstart-1, strand)
    out=[]
    with open(PANEL_G) as f:
        for ln in f:
            if not ln.strip() or ln.startswith(('#','track','browser')): continue
            c,s,e,g,sc,strd = ln.rstrip("\n").split("\t")
            s=int(s); e=int(e)
            if g not in gene2info: continue
            contig, gstart0, strand = gene2info[g]
            rs = max(0, s - gstart0)
            re_ = max(rs, e - gstart0)
            out.append((contig, rs, re_, g, ".", strand))
    with open(PANEL_P,"w") as f:
        for t in out: f.write("\t".join(map(str,t))+"\n")

sh(f"gunzip -c '{FQGZ}' | minimap2 -t 4 -ax splice -uf -k14 '{FA}' - | samtools sort -m 1G -@2 -o '{BAM}' -")
sh(f"samtools index '{BAM}'")
sh(f'stringtie -L -G "{GTF}" -o "{GTF_OUT}" "{BAM}"')
sh(f'gffread -E "{GTF_OUT}" -T -o "{GTF_TMP}"')

with open(GTF_TMP) as fin, open(BED_TX,"w") as fout:
    for ln in fin:
        if not ln or ln[0]=="#": continue
        f = ln.rstrip("\n").split("\t")
        if len(f)<9 or f[2]!="transcript": continue
        m = re.search(r'transcript_id "([^"]+)"', f[8])
        tid = m.group(1) if m else "NA"
        chrom=f[0]; start=str(int(f[3])-1); end=f[4]; strand=f[6]
        fout.write("\t".join([chrom,start,end,tid,".",strand])+"\n")

sh(f'bedtools intersect -wa -wb -a "{BED_TX}" -b "{PANEL_P}" > "{BED_INT}"')

n_all   = int(sp.check_output(["bash","-lc", f"grep -c $'\\ttranscript\\t' '{GTF_TMP}' || echo 0"]).decode().strip() or 0)
n_panel = int(sp.check_output(["bash","-lc", f"wc -l < '{BED_INT}' 2>/dev/null || echo 0"]).decode().strip() or 0)
with open(TSV_COUNTS,"w") as f:
    f.write(f"{KEEP}\t{n_all}\t{n_panel}\n")

rows=[]
if os.path.exists(BED_INT) and os.path.getsize(BED_INT)>0:
    with open(BED_INT) as fin:
        for ln in fin:
            if not ln.strip(): continue
            t=ln.rstrip("\n").split("\t")
            if len(t) < 12: continue
            rows.append((KEEP, t[9], 1))
hits = pd.DataFrame(rows, columns=["sample","gene","count"]) if rows else pd.DataFrame(columns=["sample","gene","count"])
hits.to_csv(TSV_HITS, sep="\t", header=False, index=False)

if not hits.empty:
    agg = (hits.groupby("gene", as_index=False)["count"].sum().sort_values("count", ascending=False))
    agg.to_csv(TSV_AGG, sep="\t", index=False)
    plt.figure(figsize=(8,5))
    a = agg.sort_values("count").tail(15)
    plt.barh(a["gene"], a["count"])
    plt.xlabel("Transcript overlaps (panel)")
    plt.title("Top HBOC panel genes (subsampled)")
    plt.tight_layout()
    plt.savefig(PNG_OUT, dpi=160, bbox_inches="tight")
else:
    pd.DataFrame(columns=["gene","count"]).to_csv(TSV_AGG, sep="\t", index=False)

print(open(TSV_COUNTS).read().strip())
