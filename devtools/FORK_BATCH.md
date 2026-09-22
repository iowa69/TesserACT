# Developer fork diagnostics

This Linux-only, default-off facility reuses one in-memory final graph and corrected read store for sequential **diagnostic** arms from the same executable. It does not provide a persistent checkpoint or replace standalone full-panel confirmation. The normal assembler and existing benchmark driver do not enable it.

Set `TESSERACT_DEV_FORK_BATCH` to an absolute or relative TSV manifest path. Unset it for ordinary assembly. The manifest uses exact LF-delimited `arm<TAB>flags` columns:

```text
arm	flags
baseline	-
weighted	TESSERACT_WEIGHTED_ELIGIBLE_COVERAGE=1
weighted_owned	TESSERACT_WEIGHTED_ELIGIBLE_COVERAGE=1,TESSERACT_OWNED_ANCHOR_PREFIX=1
```

In the example above, replace the displayed `\t` separators with actual tab characters. Arm IDs contain letters, digits, underscores and hyphens and begin with a letter or digit. A manifest contains 1–128 arms, is at most 1 MiB, and rejects duplicate IDs or flags. Flag values are exactly `0` or `1`; `-` means all seven optional flags are absent. A fresh output root is required: existing arm paths, batch metadata and ordinary assembly outputs are rejected. An exclusive `.fork_batch_claim` is created before preprocessing, preventing concurrent batches from sharing a root. The claim remains after success or failure; use a fresh directory for retries.

Only these flags may vary:

- `TESSERACT_ROUTE_DISTANCE`
- `TESSERACT_WEIGHTED_RESOLVER_COVERAGE`
- `TESSERACT_WEIGHTED_ELIGIBLE_COVERAGE`
- `TESSERACT_EXACT_READ_THREADS`
- `TESSERACT_OWNED_ANCHOR_PREFIX`
- `TESSERACT_EXCLUDE_SHARED_REPEAT_SUPPORT`
- `TESSERACT_OBSERVED_GRAPH_COVERAGE`

Each child independently clears all seven before applying its row, even if the parent environment set one. All other effective options and upstream environment values remain fixed. EC, mate rescue, ladder/carry, graph simplification and graph gap-closure flags cannot vary. The observed-coverage experiment is first consumed after the fork boundary; it measures the unchanged final graph from corrected reads. This live seven-flag facility is distinct from the earlier frozen six-flag diagnostic binary. No upstream construction flag is introduced by the batching facility.

Model/organism/IS options, external alignment polishing, disabled repeat resolution, unpaired input, `TESSERACT_JOIN_DUMP`, and the external model-cluster input `TESSERACT_PLASMID_CLUSTERS` are rejected. The latter is an input, not an output sink; it is excluded as an unsupported external model dependency. All inspected environment-controlled diagnostic output otherwise goes to stderr and is isolated per child.

After the final graph is built, the parent verifies exactly one OS thread and flushes stdio before each fork. One child runs at a time. A child writes normal outputs to `<root>/<arm>/`, logs to `diagnostic.log`, and provenance to `fork_child.json`. The parent writes only `fork_batch.json` and exits through an explicit batch-completion path; it emits no parent FASTA or ordinary assembly summary. Nonzero child exits/signals and cancellation fail the batch, and incomplete arms remain explicit. SIGINT/SIGTERM are recorded by a signal handler and forwarded synchronously while the child is known unreaped. A blocked parent termination signal or nondefault/nonwaitable SIGCHLD disposition is rejected. Children install default, unblocked termination signals; parent-death handling requests termination if the parent disappears. SIGKILL can leave metadata marked running; it never marks incomplete output complete.

Provenance records the real parent invocation, exact manifest/executable SHA-256, effective option values, parent `TESSERACT_*` environment, per-child flag values, shared state dimensions, statuses and PIDs. The common state identity means the unchanged in-memory parent; it is not falsely described as a serialized or cryptographically hashed read checkpoint. Source/build hashes should accompany the frozen diagnostic executable externally. Completion records normal successful output writes, exit status and expected file presence; sidecars do not independently hash/validate every output and must not be used as benchmark acceptance markers. An external byte-identity/output-hash gate remains required.

`report.json` and the child summary retain the actual parent command. Their elapsed field is explicitly **post-graph child time**, identified by the sidecar and log. Their peak memory is **child-process-only**, excluding the parent preprocessing high-water, and is labeled as such. Shared preprocessing time is recorded once in the parent metadata. Parent wait time, child CPU time and sampled parent-plus-child PSS are separate. These are not standalone assembler runtime or full-run peak-memory measurements. PSS sampling can miss brief peaks and is not a cgroup-enforced memory limit; retained allocator arenas and COW writes can increase memory. The batch consumes a resource slot for its entire lifetime.

Tests are synthetic and require explicit binaries:

```sh
python3 repo/tests/test_dev_fork_batch.py \
  --binary /absolute/path/to/diagnostic-binary \
  --baseline /absolute/path/to/frozen-baseline \
  --helper /absolute/path/to/test_dev_fork_batch \
  --out /fresh/test/output
```

`test_dev_fork_batch.cpp` includes the helper implementation to test private hashing and actual fork flow without production test hooks; link it against ordinary objects **excluding** `main.o`, `model_main.o` and `dev_fork_batch.o`. Build everything into a fresh work directory, never onto an executable used by a running campaign.

A reviewed diagnostic may identify a candidate worth testing. Acceptance still requires ordinary independent full assembly and every existing completeness, contiguity, correctness and size gate on the fixed complete panel.
