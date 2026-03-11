#!/usr/bin/env python3
import argparse, re
from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

TRPE_PATTS = [
    re.compile(r"\btrpE\b", re.I),
    re.compile(r"anthranilate synthase component\s*I\b", re.I),
    re.compile(r"aminodeoxychorismate/anthranilate synthase component\s*I\b", re.I),
]

TRPG_PATTS = [
    re.compile(r"\btrpG\b", re.I),
    re.compile(r"anthranilate synthase component\s*II\b", re.I),
    re.compile(r"aminodeoxychorismate/anthranilate synthase component\s*II\b", re.I),
]

def is_trpE(feat):
    qs = feat.qualifiers
    texts = []
    for k in ("gene", "product", "Name", "locus_tag"):
        v = qs.get(k, [])
        if isinstance(v, str):
            texts.append(v)
        else:
            texts.extend(v)
    s = " ; ".join(texts)
    return any(p.search(s) for p in TRPE_PATTS)

def is_trpG(feat):
    qs = feat.qualifiers
    texts = []
    for k in ("gene", "product", "Name", "locus_tag"):
        v = qs.get(k, [])
        if isinstance(v, str):
            texts.append(v)
        else:
            texts.extend(v)
    s = " ; ".join(texts)
    return any(p.search(s) for p in TRPG_PATTS)

def slice_record(rec, start, end, pad):
    start = max(0, start - pad)
    end   = min(len(rec.seq), end + pad)

    sub = SeqRecord(
        rec.seq[start:end],
        id=rec.id,
        name=rec.name,
        description=f"{rec.id}:{start+1}-{end}",
    )
    sub.annotations["molecule_type"] = rec.annotations.get("molecule_type", "DNA")

    new_feats = []
    for f in rec.features:
        fstart = int(f.location.start)
        fend   = int(f.location.end)
        if fend < start or fstart > end:
            continue
        new_loc = FeatureLocation(
            max(0, fstart - start),
            min(end - start, fend - start),
            strand=f.location.strand,
        )
        new_feats.append(SeqFeature(new_loc, type=f.type, qualifiers=f.qualifiers))
    sub.features = new_feats
    return sub, start, end

def process_file(path, outdir, pad, max_span=None):
    out_rows = []
    out_paths = []
    for rec in SeqIO.parse(str(path), "genbank"):
        # use CDSs only to avoid gene/CDS duplicates
        trpE_feats = [f for f in rec.features if f.type == "CDS" and is_trpE(f)]
        trpG_feats = [f for f in rec.features if f.type == "CDS" and is_trpG(f)]
        if not trpE_feats or not trpG_feats:
            continue

        best = None
        bestdist = None
        for e in trpE_feats:
            for g in trpG_feats:
                d = abs(int(e.location.start) - int(g.location.start))
                if bestdist is None or d < bestdist:
                    bestdist = d
                    best = (e, g)

        e, g = best
        raw_start = min(int(e.location.start), int(g.location.start))
        raw_end   = max(int(e.location.end),   int(g.location.end))
        span = raw_end - raw_start

        if max_span is not None and span > max_span:
            # skip if they are too far apart
            continue

        sub, start, end = slice_record(rec, raw_start, raw_end, pad)
        out_name = f"{Path(path).stem}_{rec.id}_trpE_trpG.gbk"
        out_file = outdir / out_name
        SeqIO.write(sub, out_file, "genbank")
        out_paths.append(out_file)

        ename = "/".join(e.qualifiers.get("gene", e.qualifiers.get("product", ["?"])))
        gname = "/".join(g.qualifiers.get("gene", g.qualifiers.get("product", ["?"])))
        out_rows.append((path.name, rec.id, start+1, end, span, ename, gname, out_name))

    return out_rows, out_paths

def main():
    ap = argparse.ArgumentParser(
        description="Extract regions between trpE (AS comp I) and trpG (AS comp II) from GenBank files."
    )
    ap.add_argument("--gbk", required=True,
                    help="A .gbk/.gbff file OR a directory to scan recursively")
    ap.add_argument("--out-dir", required=True,
                    help="Directory to write region .gbk files")
    ap.add_argument("--pad", type=int, default=1000,
                    help="Padding (bp) to add on both sides [default: 1000]")
    ap.add_argument("--max-span", type=int, default=None,
                    help="Maximum allowed distance (bp) between trpE and trpG CDSs "
                         "before padding; pairs with larger spans are skipped.")
    args = ap.parse_args()

    inpath = Path(args.gbk)
    outdir = Path(args.out_dir); outdir.mkdir(parents=True, exist_ok=True)

    if inpath.is_dir():
        inputs = [p for p in inpath.rglob("*")
                  if p.suffix.lower() in {".gbk", ".gb", ".gbff"}]
    else:
        inputs = [inpath]

    report_lines = [
        "source_file\trecord_id\tstart\tend\tspan_bp\ttrpE_label\ttrpG_label\tregion_file\n"
    ]
    kept = 0
    for p in inputs:
        rows, outs = process_file(p, outdir, args.pad, args.max_span)
        for r in rows:
            report_lines.append("\t".join(map(str, r)) + "\n")
        kept += len(outs)

    (outdir / "trp_island_report.tsv").write_text("".join(report_lines))
    print(f"Wrote {kept} region .gbk files to: {outdir}")
    print(f"Report: {outdir/'trp_island_report.tsv'}")

if __name__ == "__main__":
    main()
