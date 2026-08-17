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

.PHONY: all native debug asan clean install uninstall test unittest check model flagcheck probe

all: $(BIN) $(MODELBIN)

# Each of these rebuilds from scratch with different flags. Sub-makes rather
# than "clean $(BIN)" prerequisites, which make is free to run in either order
# under -j and which can therefore leave no binary at all.
native:
	@$(MAKE) --no-print-directory clean
	@$(MAKE) --no-print-directory OPT="-O3 -march=native -mtune=native" all

debug:
	@$(MAKE) --no-print-directory clean
	@$(MAKE) --no-print-directory OPT="-O0 -g3 -fno-omit-frame-pointer" all

asan:
	@$(MAKE) --no-print-directory clean
	@$(MAKE) --no-print-directory OPT="-O1 -g3 -fsanitize=address,undefined -fno-omit-frame-pointer" \
	         LDFLAGS="$(LDFLAGS) -fsanitize=address,undefined" all

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

test: $(BIN)
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

unittest: $(UNITBIN)
	@$(UNITBIN)

$(UNITBIN): $(UNITSRC) $(UNITOBJS) | $(BUILDDIR)
	$(CXX) $(CXXFLAGS) -I$(SRCDIR) $(UNITSRC) $(UNITOBJS) $(LDFLAGS) $(LDLIBS) -o $@

# Everything: unit tests then the end-to-end suite.
check: unittest test flagcheck

install: $(BIN) $(MODELBIN)
	@install -d $(DESTDIR)$(PREFIX)/bin
	@install -m 755 $(BIN) $(DESTDIR)$(PREFIX)/bin/$(BIN)
	@install -m 755 $(MODELBIN) $(DESTDIR)$(PREFIX)/bin/$(MODELBIN)
	@echo "installed $(DESTDIR)$(PREFIX)/bin/$(BIN)"
	@echo "installed $(DESTDIR)$(PREFIX)/bin/$(MODELBIN)"

uninstall:
	@rm -f $(DESTDIR)$(PREFIX)/bin/$(BIN) $(DESTDIR)$(PREFIX)/bin/$(MODELBIN)

clean:
	@rm -rf $(BUILDDIR) $(BIN) $(MODELBIN)
