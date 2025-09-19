#!/usr/bin/env bash
set -euo pipefail

# ---- config ----
BASE="${BASE:-/content/drive/MyDrive/brca-targeted-analysis/brca_rna_targeted}"
KEEP="${KEEP:-ERR13137440}"

RAW_DIR="$BASE/data/raw"
PROC_DIR="$BASE/data/processed"
REFS_DIR="$BASE/refs"

RAW_GZ="$RAW_DIR/${KEEP}.fastq.gz"
MINI_GZ="$PROC_DIR/mini_${KEEP}.fastq.gz"

FA="$REFS_DIR/GRCh38.panel28.pad50k.fa"
GTF="$REFS_DIR/gencode.v46.panel28.gtf"
BED_RAW="$REFS_DIR/panel28_genes.bed"
BED_PANELSPACE="$REFS_DIR/panel28_genes.panelspace.bed"

ENA_URL="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR131/040/ERR13137440/ERR13137440.fastq.gz"

mkdir -p "$RAW_DIR" "$PROC_DIR" "$REFS_DIR"

# ---- checks for required tools ----
for t in wget seqtk python3; do
  command -v "$t" >/dev/null 2>&1 || { echo "[MISSING] $t"; exit 1; }
done

# ---- refs  ----
test -s "$FA"  || { echo "[MISSING] $FA";  exit 1; }
test -s "$GTF" || { echo "[MISSING] $GTF"; exit 1; }
test -s "$BED_RAW" || { echo "[MISSING] $BED_RAW"; exit 1; }

# ---- download raw FASTQ if missing ----
if [ ! -s "$RAW_GZ" ]; then
  echo "[DL] $ENA_URL -> $RAW_GZ"
  wget -O "$RAW_GZ" "$ENA_URL"
else
  echo "[SKIP] raw FASTQ exists: $RAW_GZ"
fi

# ---- subsample to 30k reads into .gz ----
echo "[SUBSAMPLE] 30,000 reads -> $MINI_GZ"
seqtk sample -s100 "$RAW_GZ" 30000 | gzip -c > "$MINI_GZ"

# sanity count
reads=$(gzip -cd "$MINI_GZ" | awk 'END{print int(NR/4)}')
echo "[OK] mini FASTQ: $MINI_GZ  (${reads} reads)"

# ---- build panelspace BED from panel FASTA headers + raw BED ----
# FASTA headers look like: >GENE|chr2:47295157-47437601(+)
if [ ! -s "$BED_PANELSPACE" ]; then
  echo "[BUILD] $BED_PANELSPACE"
  python3 - "$FA" "$BED_RAW" "$BED_PANELSPACE" <<'PY'
import sys,re
fa, bedg, out = sys.argv[1], sys.argv[2], sys.argv[3]

gene2info={}
with open(fa) as f:
    for ln in f:
        if ln.startswith(">"):
            h=ln[1:].strip()
            m=re.match(r'^([^|]+)\|([^:]+):(\d+)-(\d+)\(([+-])\)', h)
            if not m: continue
            gene=m.group(1); gstart=int(m.group(3)); strand=m.group(5)
            contig=h
            gene2info[gene]=(contig, gstart-1, strand)  # 0-based start

out_rows=[]
with open(bedg) as f:
    for ln in f:
        if not ln.strip() or ln[0] in "#tb": continue
        c,s,e,g,sc,strd = ln.rstrip("\n").split("\t")
        s=int(s); e=int(e)
        if g not in gene2info: continue
        contig, gstart0, strand = gene2info[g]
        rs = max(0, s - gstart0)
        re_ = max(rs, e - gstart0)
        out_rows.append((contig, rs, re_, g, ".", strand))

with open(out,"w") as fo:
    for r in out_rows:
        fo.write("\t".join(map(str,r))+"\n")

print(f"[OK] wrote {out} with {len(out_rows)} lifted features")
PY
else
  echo "[SKIP] panelspace exists: $BED_PANELSPACE"
fi

echo "[DONE] fetch + subsample + panelspace build"
