# Phase 2, workstream emit, part B (EVAL_PLAN_P2 W1)

Branch `p2/emitb` from 1.5.0 (`4d9baaf`). Every feature is behind a `TESSERACT_P2_*` flag, default
off. With every flag unset the four assembly outputs (contigs.fasta, scaffolds.fasta, scaffolds.agp,
assembly_graph.gfa) are byte-identical to 1.5.0, no new file is written, and report.json gains only the
registered block `"p2": {"enabled": 0, "counters": {...}}` with every counter 0 (EVAL_PLAN_P2 s5.2).

| flag | feature | footprint class (s5.2) | outputs |
|---|---|---|---|
| `TESSERACT_P2_ENDS=1` | contig-end audit | report-only | `ends.tsv`, report.json `p2.ends` |
| `TESSERACT_P2_TIPS_LOWERCASE=1` | unsupported tip bases in lower case (implies the audit) | case-only | contigs.fasta tips in lower case, `ends.tsv`, `p2_edits.tsv` |
| `TESSERACT_P2_SELF_QA=1` | F2 k-mer self-QA (kqa port) | report-only | report.json `p2.self_qa`, report.html content check, `p2_self_qa.{tsv,contigs.tsv,missing_hi.tsv,absent.bed}` |
| `TESSERACT_P2_PROVENANCE=1` | F10 provenance manifest | report-only | report.json `p2.provenance` |
| `TESSERACT_P2_DETECT_REPORT=1` (+ `TESSERACT_P2_DETECT_SKETCH=<sketch>`) | F9 organism detection as a report | report-only | report.json `p2.organism_detect` (never selects or applies a model) |

Counter lines, printed on every run (zeros and `enabled=0` when off):
`[p2-ends]`, `[p2-sqa]`, `[p2-prov]`, `[p2-detect]`, `[p2-readpass]`.

## What each feature measures

**One raw read pass.** ends, self-QA and detection share one pass over the raw input files at the end of
the run (one reader thread per file, `-t` counting threads; every count is a sum, so the result does not
depend on the thread count). Nothing in it changes an assembly decision: it runs after every sequence
edit, on the records exactly as contigs.fasta receives them.

**ends.tsv** (port of the second-pass track's end audit, `backward/secondpass/proto/sp2v5.py`): one row
per end of every contigs.fasta record of at least 500 bp. Copy number (record coverage over the
length-weighted median of the records >= 5 kb), raw 31-mer end depth, non-unique tail (bases covered by a
31-mer occurring more than once in contigs.fasta, within 741 bp), pair partners from reads anchored on the
end windows (12 probes, >= 2 unique 31-mer votes, as ContigEndLinks; class unique / branching /
repeat_only / none / na), and unsupported tip bases. The tip rule is sp2's: the locus pool of an end
(reads anchored on its last 400 bp plus the mates of reads pointing out of it) is counted from the raw
reads (3' trimmed below Q10); from the end inward, positions whose 31-mer has a pool count below
max(2, 0.25 x the pool's median count over the last 300 bp) are unsupported, up to the first vouched
k-mer, at most 201. Only anchor records with >= 10 pool reads are audited (`tip_status`). The raw reads
are fetched by store index (record j of R1/R2 is store read 2j/2j+1); a run with more than one library,
or a file whose record count differs from the store's, counts the pool's corrected reads instead
(`pool_source` = store_corrected).

**Lower-case tips.** Only the bytes written to contigs.fasta change case. scaffolds.fasta, the AGP, the
GFA, every statistic and report.json `assembly` keep the upper-case records. Each lower-case run is a
record prefix or suffix of exactly `unsupported_tip_bp` bases; every edit is a `p2_edits.tsv` row
(record, feature=tips_lowercase, operation, start, end, bases; 1-based inclusive). Lower case class
`unsupported_tip`: ingested as masked bases (EVAL_PLAN_P2 s3.3 item 3).

**self_qa** (port of `backward/novel/tools/kqa.cpp` with its parameters: 31-mers hash-sampled 1 in 8,
records >= 500 bp, PhiX174 NC_001422.1 as the spike-in reference, embedded in `src/p2_phix.h`):
completeness, QV from read-absent assembly k-mers, missing sequence by depth class, spike-in k-mers, and
the frozen lost-replicon alarm of L-SQA-b (missing_bp_est_gt10x >= 1000). The summary keys of
`p2_self_qa.tsv` are kqa's `.summary.tsv` keys. With the alarm, report.html qualifies its contiguity
grade ("... -- content missing", never shown as plain ok).

**provenance**: every input file (path, real path, size, mtime, md5 of the first MiB = the armx manifest
key, SHA-256 of the whole file), the binary (/proc/self/exe; md5, SHA-256), the model and every file
named by a CLI option or a TESSERACT_* text flag (md5, SHA-256), every TESSERACT_* variable (value,
parsed value, kind), argv, the non-default options, the [defaults] line, mode, k ladder, threads, and
versions (TesserACT, git commit of the build from `build/git_commit.h`, compiler, zlib). Hashing runs on
a background thread from the start of the run.

**organism_detect**: the reads are scored against the sketch exactly as `om2::scoreReadFiles` scores
them (core markers seen >= 2 times; component test checks parity); the block lists every model's score,
the best and second, `accepted` (>= 0.85), `call`, `mixture` (second >= 0.10) and `model_applied: false`.
`--organism` handling and the om2 detection paths are untouched.

## Files (for the integrator)

New: `src/p2_emitb.{h,cpp}` (flags, driver, p2 JSON block, counter lines), `src/p2_readpass.{h,cpp}`,
`src/p2_kmer.{h,cpp}`, `src/p2_ends.{h,cpp}`, `src/p2_selfqa.{h,cpp}`, `src/p2_provenance.{h,cpp}`,
`src/p2_sha256.{h,cpp}`, `src/p2_json.h`, `src/p2_phix.h`, `tests/test_p2_emitb.cpp`, this file.

Touched:
- `src/envflags.cpp` (6 flag rows), `src/envflags.h` + `.cpp` (`env::setVariables()`).
- `src/assembler.cpp`: options and the provenance thread after `emitfix::logSummary()`; the record stage
  between `computeContigStats()` and the contigs.fasta `writeFasta` (contigs.fasta is written from the
  lower-case copy when the tips flag is on); the p2 block, HTML summary and counter lines after
  `om2::logClonalCounters()`. **Integration order:** any phase-2 feature that edits, moves or renames
  contigs.fasta records (emit-A: R3 circle closure/trim, F5 spike-in moves) must run before this stage,
  so the audit sees the records as written.
- `src/report.h` (`p2Json`, `P2Html`), `src/report_json.cpp` (the `"p2"` key before `"assembly"`),
  `src/report_html.cpp` (verdict qualification and the "Content and provenance checks" list).
- `Makefile` (`build/git_commit.h` for `p2_provenance.o`), `tests/run_tests.sh` (section 21: flags-off
  block, footprint, and the -t 1/4/6 identity test with the flags off and on).

The `"p2"` report block is one object; another workstream adds its keys next to `counters` (and its own
counters inside it) rather than a second top-level key.
