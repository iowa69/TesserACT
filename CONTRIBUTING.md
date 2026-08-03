# Contributing

Thanks for looking. A few things worth knowing before you start.

## The standard this project holds itself to

Every claim about behaviour in this repository is expected to be backed by a
measurement, and the measurement is expected to be on enough data to mean
something. Several changes in this codebase looked like clear wins on five
isolates and turned out to be neutral or harmful on forty. If you are changing
anything that affects an assembly, say what you measured it on.

That cuts both ways: a change that measures worse is a useful result and the
comment explaining *why* it was reverted is worth keeping.

## Before opening a pull request

```sh
make check      # unit tests, end-to-end tests, and the flag suite
```

All three must pass. `make asan` builds with AddressSanitizer and
UndefinedBehaviorSanitizer if you are chasing something subtle.

The three suites cover different things:

- `make unittest` — k-mer encoding at every k that spans a word boundary, the
  open-addressed table against `std::unordered_map`, the banded identity, and
  the bidirected graph invariant after every simplification operation.
- `make test` — end-to-end on synthetic genomes it generates itself, including
  that a repeat the fragments *cannot* span is left unjoined rather than
  guessed at.
- `make flagcheck` — every documented flag parses, changes behaviour, and
  rejects nonsense. A flag that is documented but silently ignored fails
  nothing at runtime, which is why this suite exists.

## Style

Match the surrounding code. Comments explain *why*, not *what* — the code
already says what it does. If a constant is load-bearing, the comment should
say what happens when it moves, ideally with the number that was measured.

## Reporting a problem

An assembly that comes out wrong is much easier to diagnose with:

- the exact command line,
- the tessera version (`tessera --version`),
- `report.json` from the output directory,
- read length, approximate coverage, and expected genome size.

If the reads can be shared, that helps enormously. If they cannot — which is
often the case for clinical isolates — `report.json` plus the command line is
usually enough to narrow it down.

## Licence

MIT. By contributing you agree your contribution is licensed under it.
