#!/usr/bin/env python3
"""Build a curated plasmid reference from PLSDB 2025 and the MOB-suite curated set.

PLSDB is comprehensive and noisy. It carries 72,556 records with a median of 53.7 kb,
and 3,908 of them (5.4%) exceed 300 kb -- a tail that includes a 5,369,247 bp
"K. pneumoniae plasmid" and two E. coli "plasmids" near 4.9 Mb. Those are chromosomes.
Left in, they put whole chromosomes of the very organisms being modelled into the
plasmid table, and the replicon classifier then has every reason to call chromosomal
contigs plasmid.

MOB-suite's curated set is the corroborating authority: 23,700 records, median 43.6 kb,
and a maximum of 399,913 bp -- nothing above 400 kb survived its curation.

Rule:
  * PLSDB record <= 300 kb                      -> keep
  * PLSDB record  > 300 kb, accession in MOB-suite -> keep (corroborated)
  * PLSDB record  > 300 kb, not corroborated    -> DROP, and record why
  * every MOB-suite record                      -> keep
Deduplicated by versionless accession, PLSDB preferred when both carry one.
"""
import sys

MAXLEN = 300_000


def read_fasta(path):
    name, seq = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq)
                name, seq = line[1:].rstrip("\n"), []
            else:
                seq.append(line.strip())
    if name is not None:
        yield name, "".join(seq)


def acc(name):
    return name.split()[0].split(".")[0]


def main(plsdb, mob, out, report):
    mob_accs = set()
    kept, dropped = {}, []
    nmob = 0
    for name, seq in read_fasta(mob):
        a = acc(name)
        mob_accs.add(a)
        kept[a] = (name, seq, "mobsuite")
        nmob += 1
    sys.stderr.write(f"MOB-suite: {nmob:,} records\n")

    npl = nkeep = ncorr = 0
    for name, seq in read_fasta(plsdb):
        npl += 1
        a = acc(name)
        if len(seq) > MAXLEN:
            if a in mob_accs:
                ncorr += 1
                kept[a] = (name, seq, "plsdb>300kb corroborated")
            else:
                dropped.append((a, len(seq), name[:100]))
                continue
        else:
            nkeep += 1
            kept[a] = (name, seq, "plsdb")
    sys.stderr.write(f"PLSDB: {npl:,} records, {nkeep:,} kept <=300kb, "
                     f"{ncorr:,} kept >300kb corroborated, {len(dropped):,} dropped\n")

    with open(out, "w") as fh:
        for a, (name, seq, src) in sorted(kept.items()):
            fh.write(f">{name.split()[0]}\n")
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i + 80] + "\n")
    with open(report, "w") as fh:
        fh.write("accession\tlength\theader\n")
        for a, ln, nm in sorted(dropped, key=lambda t: -t[1]):
            fh.write(f"{a}\t{ln}\t{nm}\n")
    sys.stderr.write(f"curated: {len(kept):,} records -> {out}\n")
    sys.stderr.write(f"dropped list -> {report}\n")


if __name__ == "__main__":
    main(*sys.argv[1:5])
