#!/usr/bin/env bash
#SBATCH -J syr_fastani_phylogroups
#SBATCH -o logs/%x.%j.out
#SBATCH -e logs/%x.%j.err
#SBATCH -p mpag-np
#SBATCH -A mpag-np
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=04:00:00

set -euo pipefail

# -----------------------------
# Inputs
# -----------------------------
IN_DIR="input_fasta_syr"
REF_GENOMES_DIR="data/genomes"

ENV_PREFIX="" # path to env with fastANI
FASTANI="${ENV_PREFIX}/bin/fastANI"

OUT_ROOT="phylogroup_fastani_out"
ANI_DIR="${OUT_ROOT}/fastani"
LOG_DIR="logs"

mkdir -p "${OUT_ROOT}" "${ANI_DIR}" "${LOG_DIR}"

# -----------------------------
# Sanity checks
# -----------------------------
[[ -d "${IN_DIR}" ]] || { echo "[ERROR] Missing IN_DIR: ${IN_DIR}" >&2; exit 1; }
[[ -d "${REF_GENOMES_DIR}" ]] || { echo "[ERROR] Missing REF_GENOMES_DIR: ${REF_GENOMES_DIR}" >&2; exit 1; }
[[ -x "${FASTANI}" ]] || { echo "[ERROR] fastANI not executable at: ${FASTANI}" >&2; exit 1; }

echo "Using fastANI: ${FASTANI}"
"${FASTANI}" --version || true

# Build refs list
REF_LIST="${ANI_DIR}/refs.fofn"
find "${REF_GENOMES_DIR}" -maxdepth 1 -type f -name "*.fna" | sort > "${REF_LIST}"
echo "[info] N reference FASTAs: $(wc -l < "${REF_LIST}")"
[[ "$(wc -l < "${REF_LIST}")" -gt 0 ]] || { echo "[ERROR] No refs found in ${REF_GENOMES_DIR}" >&2; exit 2; }

# Build query list
QUERY_LIST="${ANI_DIR}/queries.fofn"
find "${IN_DIR}" -type f \( -name "*.fa" -o -name "*.fna" -o -name "*.fasta" \) | sort > "${QUERY_LIST}"
echo "[info] N query FASTAs: $(wc -l < "${QUERY_LIST}")"
[[ "$(wc -l < "${QUERY_LIST}")" -gt 0 ]] || { echo "[ERROR] No FASTAs found in ${IN_DIR}" >&2; exit 3; }

# Run fastANI
ANI_OUT="${ANI_DIR}/fastani_query_vs_refs.tsv"
echo "[run] fastANI..."
"${FASTANI}" \
  --ql "${QUERY_LIST}" \
  --rl "${REF_LIST}" \
  -o "${ANI_OUT}" \
  --threads "${SLURM_CPUS_PER_TASK}"

echo "[done] fastANI wrote: ${ANI_OUT}"

# Best hit per query + phylogroup call (>=95% else unassigned)
BEST_HITS="${ANI_DIR}/best_hit_per_query.tsv"
ASSIGNMENTS="${OUT_ROOT}/phylogroup_assignments.tsv"

awk -F'\t' '
  BEGIN{OFS="\t"}
  {q=$1; r=$2; ani=$3+0; if (!(q in best) || ani>best[q]) {best[q]=ani; bref[q]=r}}
  END{print "query_fasta","best_ref","best_ani" > out; for (q in best) print q,bref[q],best[q] >> out}
' out="${BEST_HITS}" "${ANI_OUT}"

awk -F'\t' '
  BEGIN{OFS="\t"}
  NR==1{next}
  {
    q=$1; r=$2; ani=$3+0;

    # genome_id from query basename without extension
    qb=q; sub(/^.*\//,"",qb);
    gid=qb; sub(/\.fasta$/,"",gid); sub(/\.fna$/,"",gid); sub(/\.fa$/,"",gid);

    # phylogroup from ref basename: <pg>__GCF_... .fna
    rb=r; sub(/^.*\//,"",rb);
    split(rb,a,"__"); pg=a[1];

    call=(ani>=95 ? pg : "unassigned");
    print gid,q,ani,pg,call
  }
' "${BEST_HITS}" | sort -t$'\t' -k1,1 > "${ASSIGNMENTS}"

tmp="${ASSIGNMENTS}.tmp"
{
  echo -e "genome_id\tquery_fasta\tbest_ani\tbest_ref_phylogroup\tassigned_phylogroup_95pct"
  cat "${ASSIGNMENTS}"
} > "${tmp}"
mv "${tmp}" "${ASSIGNMENTS}"

echo "[done] Assignments: ${ASSIGNMENTS}"