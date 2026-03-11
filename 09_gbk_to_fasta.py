#!/usr/bin/env python3
'''
USAGE:

python 09_gbk_to_fasta.py \
    --gbk-dir trp_islands_gbk/ \
    --out-dir trp_islands_fasta/

'''

import argparse
from pathlib import Path
from Bio import SeqIO

def gbk_to_fasta(gbk_path, outdir):
    recs = list(SeqIO.parse(str(gbk_path), "genbank"))
    if not recs:
        return None

    # use the single record (your trp-island extractor produces 1 per file)
    rec = recs[0]

    outpath = outdir / (gbk_path.stem + ".fasta")
    SeqIO.write(rec, outpath, "fasta")
    return outpath

def main():
    ap = argparse.ArgumentParser(
        description="Extract FASTA sequences from trp-island GenBank files."
    )
    ap.add_argument("--gbk-dir", required=True,
                    help="Directory containing trp-island .gbk files")
    ap.add_argument("--out-dir", required=True,
                    help="Directory to write .fasta files")
    args = ap.parse_args()

    gbk_dir = Path(args.gbk_dir)
    outdir = Path(args.out_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    for gbk in gbk_dir.iterdir():
        if gbk.suffix.lower() in {".gbk", ".gb"}:
            outpath = gbk_to_fasta(gbk, outdir)
            if outpath:
                print(f"Wrote: {outpath}")

if __name__ == "__main__":
    main()

