blastp \
  -query ~/blastdb/data/syr_tf.faa \
  -db ~/blastdb/data/ref10_db \
  -evalue 1e-5 \
  -max_target_seqs 6 \
  -outfmt "6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore qcovs" \
  -out ~/blastdb/data/syr_vs_ref10.tsv \
  2> ~/blastdb/data/syr_vs_ref10.blastp.err