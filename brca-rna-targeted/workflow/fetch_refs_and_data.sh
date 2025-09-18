#!/usr/bin/env bash
# fetch_refs_and_data.sh
# Download refs (GRCh38 + GENCODE v46), build 28-gene panel BED,
# fetch one ENA FASTQ for $KEEP, subsample to <1 GB, and compress to .xz (python lzma).
set -euo pipefail

# --------- Params (override via env) ----------
BASE="${BASE:-$PWD/brca-rna-targeted}"
KEEP="${KEEP:-ERR13137440}"          # ENA run accession
NREADS="${NREADS:-300000}"           # subsample reads
RUN_URL="${RUN_URL:-}"               # optional: direct FASTQ URL; else resolve via ENA API
INSTALL="${INSTALL:-0}"              # set 1 to apt-get curl/wget/seqtk/pigz if needed

# --------- Tools check / optional install -----
need() { command -v "$1" >/dev/null 2>&1; }
if [ "$INSTALL" = "1" ]; then
  if need apt-get; then
    sudo apt-get -qq update
    sudo apt-get -qq install -y curl wget seqtk pigz || true
  fi
fi
for t in curl wget seqtk python3; do
  need "$t" || { echo "[ERR] missing '$t' (set INSTALL=1 to attempt apt-get)"; exit 1; }
done

# --------- Layout -----------------------------
RAW="$BASE/data/raw"
PROC="$BASE/data/processed"
REF="$BASE/data/refs"
ALIGN="$BASE/align"
ISO="$BASE/isoforms"
RES="$BASE/results"
mkdir -p "$RAW" "$PROC" "$REF" "$ALIGN" "$ISO" "$RES" "$BASE/workflow/outputs"

# --------- Refs (GRCh38 + GENCODE v46) --------
FA_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/GRCh38.primary_assembly.genome.fa.gz"
GTF_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.annotation.gtf.gz"
FA="$REF/GRCh38.primary_assembly.genome.fa"
GTF="$REF/gencode.v46.annotation.gtf"

if [ ! -s "$FA" ]; then
  echo "[REF] downloading GRCh38 fasta"
  wget -q -O "$FA.gz" "$FA_URL"
  pigz -d -f "$FA.gz"
fi
if [ ! -s "$GTF" ]; then
  echo "[REF] downloading GENCODE v46 gtf"
  wget -q -O "$GTF.gz" "$GTF_URL"
  pigz -d -f "$GTF.gz"
fi
echo "[REF] done: $(basename "$FA"), $(basename "$GTF")"

# --------- 28-gene panel BED (robust, via Python) ----
PANEL_TXT="$REF/panel28_genes.txt"
cat > "$PANEL_TXT" <<'TXT'
BRCA1
BRCA2
ATM
ATR
BARD1
BRIP1
CDH1
CHEK2
MRE11
NBN
PALB2
PTEN
RAD50
RAD51
RAD51C
RAD51D
RECQL
TP53
XRCC2
XRCC3
FANCA
FANCC
FANCM
MSH2
MSH6
MLH1
PMS2
EPCAM
TXT

PANEL_BED_RAW="$REF/panel28_genes.raw.bed"
PANEL_BED="$REF/panel28_genes.bed"
python3 - "$GTF" "$PANEL_TXT" "$PANEL_BED_RAW" "$PANEL_BED" <<'PY'
import sys,re,shutil
gtf,panel_txt,raw_bed,bed = sys.argv[1:]
genes={x.strip() for x in open(panel_txt) if x.strip()}
n=0
with open(gtf) as fin, open(raw_bed,"w") as fout:
  for ln in fin:
    if not ln or ln[0]=="#": continue
    f=ln.rstrip("\n").split("\t")
    if len(f)<9 or f[2]!="gene": continue
    m=re.search(r'gene_name "([^"]+)"', f[8])
    if not m: continue
    name=m.group(1)
    if name in genes:
      chrom=f[0]; start=str(int(f[3])-1); end=f[4]; strand=f[6]
      fout.write("\t".join([chrom,start,end,name,".",strand])+"\n")
      n+=1
shutil.copy(raw_bed, bed)
print(f"[PANEL] rows:", n)
PY
echo "[PANEL] wrote: $PANEL_BED"

# --------- Resolve FASTQ URL (if not provided) ----
if [ -z "$RUN_URL" ]; then
  echo "[ENA] resolving FASTQ URL for $KEEP"
  # ENA portal (tsv; field 'fastq_ftp'); pick first URL if multiple
  FASTQ_FTP=$(curl -s "https://www.ebi.ac.uk/ena/portal/api/search?result=read_run&query=run_accession=${KEEP}&fields=fastq_ftp&format=tsv&limit=0" \
            | awk 'NR==2{print $1}' | tr ';' '\n' | head -n1)
  if [ -z "$FASTQ_FTP" ]; then
    echo "[ERR] could not resolve fastq_ftp for $KEEP"; exit 1
  fi
  RUN_URL="ftp://${FASTQ_FTP}"
fi
echo "[ENA] using FASTQ URL: $RUN_URL"

# --------- Download raw FASTQ ------------------------
RAW_FASTQ="$RAW/${KEEP}.fastq.gz"
if [ ! -s "$RAW_FASTQ" ]; then
  echo "[DL] $RUN_URL -> $RAW_FASTQ"
  wget -q -O "$RAW_FASTQ" "$RUN_URL"
else
  echo "[SKIP] raw exists: $RAW_FASTQ"
fi

# --------- Subsample to processed mini FASTQ (.gz) ----
MINI_GZ="$PROC/mini_${KEEP}.fastq.gz"
if [ ! -s "$MINI_GZ" ]; then
  echo "[SUBSAMPLE] $NREADS reads -> $MINI_GZ"
  seqtk sample -s100 "$RAW_FASTQ" "$NREADS" | pigz -c > "$MINI_GZ"
else
  echo "[SKIP] processed exists: $MINI_GZ"
fi

# --------- Convert .gz -> .xz using Python lzma ----------
MINI_XZ="$PROC/mini_${KEEP}.fastq.xz"
if [ ! -s "$MINI_XZ" ]; then
  echo "[XZ] converting $MINI_GZ -> $MINI_XZ (python lzma)"
  python3 - "$MINI_GZ" "$MINI_XZ" <<'PY'
import sys,gzip,lzma,shutil,os
src,dst = sys.argv[1], sys.argv[2]
with gzip.open(src,"rb") as fin, lzma.open(dst,"wb", preset=6) as fout:
  shutil.copyfileobj(fin,fout)
print("[XZ] wrote:", dst, "size_MB=", round(os.path.getsize(dst)/1024/1024,1))
PY
else
  echo "[SKIP] xz exists: $MINI_XZ"
fi

# --------- Summary --------------------------------------
echo
echo "[SUMMARY]"
du -h "$REF" | tail -n1
du -h "$RAW" | tail -n1
du -h "$PROC" | tail -n1
ls -lh "$PROC" | sed -n '1,999p'
echo "[DONE]"
