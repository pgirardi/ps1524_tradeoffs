#!/usr/bin/env bash
set -euo pipefail

# Output folder (will be created)
OUT_DIR="data"
GENOMES_DIR="${OUT_DIR}/genomes"
mkdir -p "${GENOMES_DIR}"

# Table 2 mapping: accession -> phylogroup
cat > "${OUT_DIR}/phylogroup_refs.tsv" << 'EOF'
accession	phylogroup
GCF_000172895.1	1a
GCF_001910465.1	1b
GCF_000145825.2	2a
GCF_003698965.1	2b
GCF_000177515.1	2c
GCF_003205905.1	2d
GCF_000012205.1	3
GCF_000156995.2	4
GCF_016599635.1	5
GCF_008692855.1	6
GCF_000452485.1	7
GCF_000452825.1	9
GCF_000452665.1	10a
GCF_000452785.1	10b
GCF_900104015.1	11
GCF_000452865.1	13
EOF

# Tool checks
command -v datasets >/dev/null 2>&1 || { echo "[ERROR] datasets not found in PATH"; exit 1; }
command -v unzip >/dev/null 2>&1 || { echo "[ERROR] unzip not found"; exit 1; }

echo "datasets: $(datasets --version || true)"

# Download each accession and extract *_genomic.fna
tail -n +2 "${OUT_DIR}/table2_phylogroup_refs.tsv" | while IFS=$'\t' read -r acc pg; do
  out_fna="${GENOMES_DIR}/${pg}__${acc}.fna"
  if [[ -s "${out_fna}" ]]; then
    echo "[skip] ${out_fna}"
    continue
  fi

  tmp="$(mktemp -d)"
  zip_path="${tmp}/${acc}.zip"

  echo "[get] ${acc} (phylogroup ${pg})"
  datasets download genome accession "${acc}" --include genome --filename "${zip_path}"

  unzip -q "${zip_path}" -d "${tmp}/unz"
  fna="$(find "${tmp}/unz" -type f -name "*_genomic.fna" -print -quit)"

  if [[ -z "${fna}" || ! -s "${fna}" ]]; then
    echo "[ERROR] Could not find *_genomic.fna inside datasets zip for ${acc}" >&2
    rm -rf "${tmp}"
    exit 2
  fi

  cp "${fna}" "${out_fna}"
  rm -rf "${tmp}"
done

echo
echo "[done] Reference genomes saved to: ${GENOMES_DIR}"
echo "List:"
ls -lh "${GENOMES_DIR}"
