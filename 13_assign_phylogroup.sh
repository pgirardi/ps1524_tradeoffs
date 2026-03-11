# based on bit score
sort -t $'\t' -k1,1 -k12,12nr ~/blastdb/data/syr_vs_ref10.tsv \
| awk -F $'\t' '!seen[$1]++' \
> ~/blastdb/data/syr_vs_ref10.best.tsv