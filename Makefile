# TesserACT - de novo short-read assembler

BIN       := tesseract-asm
MODELBIN  := tesseract-model
SRCDIR    := src
BUILDDIR  := build
PREFIX    ?= /usr/local

SOURCES   := $(wildcard $(SRCDIR)/*.cpp)
ALLOBJS   := $(patsubst $(SRCDIR)/%.cpp,$(BUILDDIR)/%.o,$(SOURCES))
# model_main.o owns a second main(); it links only into the model builder.
OBJECTS   := $(filter-out $(BUILDDIR)/model_main.o,$(ALLOBJS))
MODELOBJS := $(filter-out $(BUILDDIR)/main.o,$(ALLOBJS))
DEPS      := $(ALLOBJS:.o=.d)

# The unit tests link against everything but the two files that own main().
UNITSRC   := tests/test_units.cpp
UNITBIN   := $(BUILDDIR)/test_units
UNITOBJS  := $(filter-out $(BUILDDIR)/main.o $(BUILDDIR)/model_main.o,$(ALLOBJS))

CXX       ?= g++
CXXSTD    := -std=c++17
WARN      := -Wall -Wextra -Wno-unused-parameter
OPT       ?= -O3
CXXFLAGS  += $(CXXSTD) $(WARN) $(OPT) -pthread -MMD -MP
LDFLAGS   += -pthread
LDLIBS    += -lz

# Development-only driver for the join stage. Lives outside src/ so the wildcard
# above never sees its main(), and is never built by `all`.
PROBESRC  := devtools/join_probe.cpp
PROBEBIN  := $(BUILDDIR)/join_probe
PROBEOBJS := $(filter-out $(BUILDDIR)/main.o $(BUILDDIR)/model_main.o,$(ALLOBJS))

.PHONY: componenttest all native debug asan clean install uninstall test unittest check model flagcheck probe \
        envguard pytest recipecheck

# The model builder is deliberately NOT part of `all` or `install`. A model's value
# is in how its panel was assembled and what was withheld from it; a file of the
# right shape built from an arbitrary panel yields confident joins with nothing
# behind them. Users get bundled, checksummed models via tesseract-get-models and
# select them with --organism. Build the tool with `make model` if you are the one
# curating the panels.
all: $(BIN)

# Each of these rebuilds from scratch with different flags. Sub-makes rather
# than "clean $(BIN)" prerequisites, which make is free to run in either order
# under -j and which can therefore leave no binary at all.
native:
	@$(MAKE) --no-print-directory clean
	@$(MAKE) --no-print-directory OPT="-O3 -march=native -mtune=native" all

debug:
	@$(MAKE) --no-print-directory clean
	@$(MAKE) --no-print-directory OPT="-O0 -g3 -fno-omit-frame-pointer" all

# The model builder is linked here too: `make test` needs it (run_tests.sh builds a track
# model), and linking it later from these objects without the sanitizer flags fails.
# -fno-sanitize-recover: an UndefinedBehaviorSanitizer report ends the run. Recoverable, it
# went to stderr, which the end-to-end suite captures and discards, and the run exited 0.
asan:
	@$(MAKE) --no-print-directory clean
	@$(MAKE) --no-print-directory OPT="-O1 -g3 -fsanitize=address,undefined -fno-sanitize-recover=undefined -fno-omit-frame-pointer" \
	         LDFLAGS="$(LDFLAGS) -fsanitize=address,undefined" all model

$(BIN): $(OBJECTS)
	$(CXX) $(OBJECTS) $(LDFLAGS) $(LDLIBS) -o $@

model: $(MODELBIN)

$(MODELBIN): $(MODELOBJS)
	$(CXX) $(MODELOBJS) $(LDFLAGS) $(LDLIBS) -o $@

$(BUILDDIR)/%.o: $(SRCDIR)/%.cpp | $(BUILDDIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BUILDDIR):
	@mkdir -p $(BUILDDIR)

-include $(DEPS)
-include $(BUILDDIR)/test_units.d

# run_tests.sh builds a track model with tesseract-model (test 15c fails without it), and
# `check` reaches `test` before `flagcheck`, so the model builder is a prerequisite here.
test: $(BIN) $(MODELBIN)
	@bash tests/run_tests.sh

flagcheck: $(BIN) $(MODELBIN)
	@bash tests/check_flags.sh ./$(BIN) ./$(MODELBIN)

PIDXSRC := devtools/plasmid_index.cpp
PIDXBIN := $(BUILDDIR)/plasmid_index

probe: $(PROBEBIN) $(PIDXBIN)

$(PIDXBIN): $(PIDXSRC) $(PROBEOBJS) | $(BUILDDIR)
	$(CXX) $(CXXFLAGS) -I$(SRCDIR) $(PIDXSRC) $(PROBEOBJS) $(LDFLAGS) $(LDLIBS) -o $@

$(PROBEBIN): $(PROBESRC) $(PROBEOBJS) | $(BUILDDIR)
	$(CXX) $(CXXFLAGS) -I$(SRCDIR) $(PROBESRC) $(PROBEOBJS) $(LDFLAGS) $(LDLIBS) -o $@

# Tests run with every ambient TESSERACT_* removed, so an exported experiment environment can
# neither fail them nor mask what they check. Each test main() also clears the environment
# itself (tests/test_env.h); `envguard` below checks that the binaries do so on their own.
TESTENV = env $$(env | sed -n 's/^\(TESSERACT_[A-Za-z0-9_]*\)=.*/-u \1/p')

unittest: $(UNITBIN)
	@$(TESTENV) $(UNITBIN)

# Standalone component tests. Each owns main() and links against the library objects.
# Until 1.3.0 none of these had a build target, so none of them ever ran in `make check`.
# test_dev_fork_batch is a helper driven by tests/test_dev_fork_batch.py, not a test.
COMPSRC   := $(filter-out tests/test_units.cpp tests/test_dev_fork_batch.cpp,$(wildcard tests/test_*.cpp))
COMPBINS  := $(patsubst tests/%.cpp,$(BUILDDIR)/%,$(COMPSRC))

$(BUILDDIR)/test_%: tests/test_%.cpp $(UNITOBJS) | $(BUILDDIR)
	$(CXX) $(CXXFLAGS) -I$(SRCDIR) $< $(UNITOBJS) $(LDFLAGS) $(LDLIBS) -o $@

componenttest: $(COMPBINS)
	@set -e; for t in $(COMPBINS); do $(TESTENV) $$t; done

# Hostile-environment guard: every unit and component test, run directly (no TESTENV), must
# still pass (a) under the K2 experiment arm's environment, which aborted three resolver tests
# before they cleared it, and (b) with every TESSERACT_* name the sources read set to 1.
SRCFLAGS  = $(shell grep -ohE '"TESSERACT_[A-Z0-9_]+"' $(SRCDIR)/*.cpp $(SRCDIR)/*.h | tr -d '"' | sort -u)
K2ENV     := TESSERACT_COMMON_PREFIX=0 TESSERACT_MIN_FALLBACK_DEST=1000000000 \
             TESSERACT_REQUIRE_SUPPORT_SINGLE=1 TESSERACT_EXCLUDE_SHARED_REPEAT_SUPPORT=1
envguard: $(UNITBIN) $(COMPBINS)
	@set -e; mkdir -p $(BUILDDIR)/envguard; \
	for e in k2 all1; do \
	  if [ $$e = k2 ]; then vars="$(K2ENV)"; else vars="$(addsuffix =1,$(SRCFLAGS))"; fi; \
	  for t in $(UNITBIN) $(COMPBINS); do \
	    log=$(BUILDDIR)/envguard/$$(basename $$t).$$e.log; \
	    env $$vars $$t > $$log 2>&1 || { echo "envguard: $$t failed under the $$e environment:"; tail -5 $$log; exit 1; }; \
	  done; \
	done; \
	echo "envguard: $(words $(UNITBIN) $(COMPBINS)) tests pass under the K2 environment and with $(words $(SRCFLAGS)) source flags set to 1"

# The python integration tests. The fork-batch helper #includes dev_fork_batch.cpp, so it
# links without dev_fork_batch.o (the test_% pattern rule would give duplicate symbols).
# PYTEST_BASELINE defaults to the binary under test, and then the A/B identity checks prove
# self-consistency only; point it at a release binary for a real A/B comparison.
FORKOBJS  := $(filter-out $(BUILDDIR)/dev_fork_batch.o,$(UNITOBJS))
FORKBIN   := $(BUILDDIR)/test_dev_fork_batch
PYTEST_BASELINE ?= $(CURDIR)/$(BIN)
$(FORKBIN): tests/test_dev_fork_batch.cpp $(FORKOBJS) | $(BUILDDIR)
	$(CXX) $(CXXFLAGS) -I$(SRCDIR) $< $(FORKOBJS) $(LDFLAGS) $(LDLIBS) -o $@

pytest: $(BIN) $(FORKBIN)
	@d=$$(mktemp -d "$${TMPDIR:-/tmp}/tesseract-pytest.XXXXXXXX"); set -e; \
	python3 tests/test_ec_memory_budget.py $(CURDIR)/$(BIN) $(PYTEST_BASELINE) $$d/ec > $$d/ec.log; \
	python3 tests/test_dev_fork_batch.py --binary $(CURDIR)/$(BIN) --baseline $(PYTEST_BASELINE) \
	        --helper $(CURDIR)/$(FORKBIN) --out $$d/fb; \
	rm -rf $$d; echo "python integration tests passed (EC memory budget, dev fork batch)"

# The conda recipe, checked against this tree: version, installed files, test commands and
# the package test, with the compile replaced by the binary already built here.
recipecheck: $(BIN)
	@bash tests/check_conda_recipe.sh

# Header dependencies of the test binaries (tests/test_env.h and the library headers).
-include $(patsubst %,%.d,$(COMPBINS) $(FORKBIN))

$(UNITBIN): $(UNITSRC) $(UNITOBJS) | $(BUILDDIR)
	$(CXX) $(CXXFLAGS) -I$(SRCDIR) $(UNITSRC) $(UNITOBJS) $(LDFLAGS) $(LDLIBS) -o $@

# Everything: unit tests then the end-to-end suite.
check: unittest componenttest envguard test flagcheck pytest recipecheck

install: $(BIN)
	@install -d $(DESTDIR)$(PREFIX)/bin
	@install -m 755 $(BIN) $(DESTDIR)$(PREFIX)/bin/$(BIN)
	@echo "installed $(DESTDIR)$(PREFIX)/bin/$(BIN)"

uninstall:
	@rm -f $(DESTDIR)$(PREFIX)/bin/$(BIN) $(DESTDIR)$(PREFIX)/bin/$(MODELBIN)

clean:
	@rm -rf $(BUILDDIR) $(BIN) $(MODELBIN)
