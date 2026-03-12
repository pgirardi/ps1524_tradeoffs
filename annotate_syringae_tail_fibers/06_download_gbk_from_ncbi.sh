# make assembly accession list
grep -Rho "^>.*assembly=GCF_[0-9]\+\.[0-9]\+" input_full_genomes \
  | sed 's/.*assembly=//' \
  | sort -u \
  > assembly_accessions.txt

# download gbk files from NCBI by accession
datasets download genome accession \
  --inputfile assembly_accessions.txt \
  --include gbff

unzip ncbi_dataset.zip

# rename to match the .fna files
mkdir -p rename_work

for f in input_full_genomes/*.fna; do
  acc=$(grep -m1 "^>" "$f" | tr -d '\r' | grep -oE 'GC[AF]_[0-9]+\.[0-9]+' | head -n 1)
  if [ -n "$acc" ]; then
    printf "%s\t%s\n" "$acc" "$f"
  else
    printf "NO_ASSEMBLY\t%s\n" "$f" >&2
  fi
done | sort -u > rename_work/assembly_to_fna.tsv

# flatten/rename
mkdir -p renamed_gbk

while IFS=$'\t' read -r acc fna; do
  base=$(basename "$fna" .fna)

  gbk="ncbi_dataset/data/$acc/genomic.gbff"
  if [ ! -f "$gbk" ]; then
    gbk=$(find "ncbi_dataset/data/$acc" -type f -name "*genomic*.gbff" -print -quit 2>/dev/null)
  fi

  if [ -z "$gbk" ] || [ ! -f "$gbk" ]; then
    echo "Missing GBK for $acc ($base)" >&2
    continue
  fi

  cp "$gbk" "renamed_gbk/${base}.gbk"
done < rename_work/assembly_to_fna.tsv

