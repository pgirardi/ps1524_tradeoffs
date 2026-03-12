# tailocin_tradeoffs

pt 1: syringae tail fibers
0) Download draft and complete P. syringae from pseudomonas.com (full genomes)
1) Run seqkit stats to generate assembly_stats.tsv
2) Run qc notebook -- qc based on number of contigs, summed length of contigs
3) Compute N50
4) Run the N50 notebook -- drop anything below 100kb
5) Run dRep for 99% ANI & selecting representative genomes per cluster
6) download genbank files for subset genomes from NCBI
7) extract trpE/trpG islands
8) QC trpE/trpG islands
9) extract nucleotide fasta from GBKs
10) run pharokka/phold on fasta using snakemake. manually curate tail fibers > syringae_tail_fibers.xlsx (inferred via synteny)

pt 2: assign syringae genomes to phylogroup by ani
1) download reference strains for phylogroups
2) assign syringae genomes to phylogroup by ANI
