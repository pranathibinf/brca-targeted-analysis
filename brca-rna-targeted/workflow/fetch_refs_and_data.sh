#!/usr/bin/env bash
# fetch_refs_and_data.sh
# Download refs (GRCh38 + GENCODE v46), build 28-gene panel BED,
# fetch one ENA FASTQ for $KEEP, subsample (default 50k reads), write .xz.

set -euo pipefail

# ------------ Params ------------
BASE="${BASE:-$PWD/brca-rna-targeted}"
KEEP="${KEEP:-ERR13137440}"          # ENA run accession
NREADS="${NREADS:-50000}"            # smaller
RUN_URL="${RUN_URL:-}"               # provide to skip ENA lookup
INSTALL="${INSTALL:-0}"              # set 1 to apt-get curl/wget/seqtk if needed
XZ_PRESET="${XZ_PRESET:-3}"          # 3 = fast; 6–9 = smaller but slower

# ------------ check tools/ optional install -------
need() { command -v "$1" >/dev/null 2>&1; }
if [ "$INSTALL" = "1" ] && need apt-get; then
  sudo apt-get -qq update
  sudo apt-get -qq install -y curl wget seqtk || true
fi
for t in curl wget python3; do
  need "$t" || { echo "[ERR] missing '$t' (set INSTALL=1 to attempt apt-get)"; exit 1; }
done

# ------------ Layout ------------------------------
RAW="$BASE/data/raw"
PROC="$BASE/data/processed"
REF="$BASE/refs"                      # <— use refs/ (not data/refs)
ALIGN="$BASE/align"
ISO="$BASE/isoforms"
RES="$BASE/results"
mkdir -p "$RAW" "$PROC" "$REF" "$ALIGN" "$ISO" "$RES" "$BASE/workflow/outputs"

# ------------ Refs (GRCh38 + GENCODE v46) --------
FA_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/GRCh38.primary_assembly.genome.fa.gz"
GTF_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.annotation.gtf.gz"
FA="$REF/GRCh38.primary_assembly.genome.fa"
GTF="$REF/gencode.v46.annotation.gtf"

if [ ! -s "$FA" ]; then
  echo "[REF] downloading GRCh38 fasta"
  wget -q -O "$FA.gz" "$FA_URL"
  python3 - <<PY
import gzip,shutil
with gzip.open("$FA.gz","rb") as fin, open("$FA","wb") as fout:
    shutil.copyfileobj(fin,fout)
PY
fi

if [ ! -s "$GTF" ]; then
  echo "[REF] downloading GENCODE v46 gtf"
  wget -q -O "$GTF.gz" "$GTF_URL"
  python3 - <<PY
import gzip,shutil
with gzip.open("$GTF.gz","rb") as fin, open("$GTF","wb") as fout:
    shutil.copyfileobj(fin,fout)
PY
fi
echo "[REF] ready: $(basename "$FA"), $(basename "$GTF")"

# ------------ 28-gene panel BED  --
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

# ------------ Resolve FASTQ URL (if not provided) -
if [ -z "$RUN_URL" ]; then
  echo "[ENA] resolving FASTQ URL for $KEEP"
  FASTQ_FTP=$(curl -s "https://www.ebi.ac.uk/ena/portal/api/search?result=read_run&query=run_accession=${KEEP}&fields=fastq_ftp&format=tsv&limit=0" \
            | awk 'NR==2{print $1}' | tr ';' '\n' | head -n1)
  [ -n "$FASTQ_FTP" ] || { echo "[ERR] could not resolve fastq_ftp for $KEEP"; exit 1; }
  RUN_URL="ftp://${FASTQ_FTP}"
fi
echo "[ENA] using FASTQ URL: $RUN_URL"

# ------------ Download raw FASTQ ------------------
RAW_FASTQ="$RAW/${KEEP}.fastq.gz"
if [ ! -s "$RAW_FASTQ" ]; then
  echo "[DL] $RUN_URL -> $RAW_FASTQ"
  wget -q -O "$RAW_FASTQ" "$RUN_URL"
else
  echo "[SKIP] raw exists: $RAW_FASTQ"
fi

# ------------ Subsample directly to .xz -----------
MINI_XZ="$PROC/mini_${KEEP}.fastq.xz"
if [ ! -s "$MINI_XZ" ]; then
  echo "[SUBSET] ${NREADS} reads -> $MINI_XZ (xz preset $XZ_PRESET)"
  python3 - "$RAW_FASTQ" "$MINI_XZ" "$NREADS" "$XZ_PRESET" <<'PY'
import sys,gzip,lzma,os
src,dst,nreads,preset = sys.argv[1],sys.argv[2],int(sys.argv[3]),int(sys.argv[4])
reads=0
with gzip.open(src,"rt",encoding="utf-8",errors="ignore") as fin, \
     lzma.open(dst,"wt",encoding="utf-8",preset=preset) as fout:
  while reads < nreads:
    h=fin.readline()
    if not h: break
    s=fin.readline(); plus=fin.readline(); q=fin.readline()
    if not q: break
    if not h.startswith("@"): continue
    fout.write(h); fout.write(s); fout.write(plus); fout.write(q)
    reads += 1
size = os.path.getsize(dst)/1024/1024
print(f"[WROTE] {dst} reads={reads} size_MB={size:.1f}")
PY
else
  echo "[SKIP] mini exists: $MINI_XZ"
fi

# ------------ Summary -----------------------------
echo
echo "[SUMMARY]"
du -h "$REF" | tail -n1 || true
du -h "$RAW" | tail -n1 || true
du -h "$PROC" | tail -n1 || true
ls -lh "$PROC" || true
echo "[DONE]"
