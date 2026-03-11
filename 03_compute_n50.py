from pathlib import Path

def n50(lengths):
    total = sum(lengths)
    half = total / 2
    running = 0
    for L in sorted(lengths, reverse=True):
        running += L
        if running >= half:
            return L
    return 0

print("file\tcontigs\tsum_len\tN50\tmax_contig")

for fasta in sorted(Path("fna").glob("*.fna")):
    lengths = []
    seq_len = 0

    with open(fasta) as fh:
        for line in fh:
            if line.startswith(">"):
                if seq_len:
                    lengths.append(seq_len)
                seq_len = 0
            else:
                seq_len += len(line.strip())
        if seq_len:
            lengths.append(seq_len)

    if not lengths:
        print(f"{fasta}\t0\t0\t0\t0")
        continue

    print(f"{fasta}\t{len(lengths)}\t{sum(lengths)}\t{n50(lengths)}\t{max(lengths)}")