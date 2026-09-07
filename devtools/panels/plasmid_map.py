#!/usr/bin/env python3
"""Map every plasmid record back to the assembly it came from.

--exclude-plasmids matches the exact record name, and a list of bare accessions
excludes nothing while still reporting success (src/model_main.cpp, and the
warning at the top of docs/MODELS_INTEGRATED.md). Withholding a held-out clone's
plasmids therefore needs this mapping: assembly -> the record names it
contributed.
"""
import glob, os, sys
raw, out = sys.argv[1], sys.argv[2]
n = 0
with open(out, "w") as fh:
    fh.write("safe_acc\taccession\tplasmid_record\tlength\n")
    for d in sorted(glob.glob(os.path.join(raw, "batch.*.d", "ncbi_dataset", "data", "GC*"))):
        acc = os.path.basename(d)
        safe = acc.replace("_", "").replace(".", "v")
        for fn in glob.glob(os.path.join(d, "*.fna")):
            name, ln = None, 0
            with open(fn) as f:
                for line in f:
                    if line.startswith(">"):
                        if name and "plasmid" in name.lower():
                            fh.write(f"{safe}\t{acc}\t{name.split()[0]}\t{ln}\n"); n += 1
                        name, ln = line[1:].rstrip("\n"), 0
                    else:
                        ln += len(line.strip())
            if name and "plasmid" in name.lower():
                fh.write(f"{safe}\t{acc}\t{name.split()[0]}\t{ln}\n"); n += 1
print(f"{n} plasmid records mapped -> {out}")
