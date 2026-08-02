# tessera - de novo short-read assembler

BIN       := tessera
MODELBIN  := tessera-model
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

.PHONY: all native debug asan clean install uninstall test unittest check model flagcheck

all: $(BIN) $(MODELBIN)

native: OPT := -O3 -march=native -mtune=native
native: clean $(BIN)

debug: OPT := -O0 -g3 -fno-omit-frame-pointer
debug: clean $(BIN)

asan: OPT := -O1 -g3 -fsanitize=address,undefined -fno-omit-frame-pointer
asan: LDFLAGS += -fsanitize=address,undefined
asan: clean $(BIN)

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

unittest: $(UNITBIN)
	@$(UNITBIN)

$(UNITBIN): $(UNITSRC) $(UNITOBJS) | $(BUILDDIR)
	$(CXX) $(CXXFLAGS) -I$(SRCDIR) $(UNITSRC) $(UNITOBJS) $(LDFLAGS) $(LDLIBS) -o $@

# Everything: unit tests then the end-to-end suite.
check: unittest test flagcheck

install: $(BIN)
	@install -d $(DESTDIR)$(PREFIX)/bin
	@install -m 755 $(BIN) $(DESTDIR)$(PREFIX)/bin/$(BIN)
	@echo "installed $(DESTDIR)$(PREFIX)/bin/$(BIN)"

uninstall:
	@rm -f $(DESTDIR)$(PREFIX)/bin/$(BIN)

clean:
	@rm -rf $(BUILDDIR) $(BIN)
