// Strict, validated reading of the TESSERACT_* environment flags.
//
// Until build_v3 every flag was read with getenv + atoi/atof at its point of use, most of
// them into a function-local static. Three things went wrong silently:
//
//   * a malformed value became 0: TESSERACT_MIN_FALLBACK_DEST=1e9 parsed as 1 (guard off),
//     TESSERACT_JOIN_TRACE=true as 0 (no trace), TESSERACT_RC_DOVETAIL=true turned a
//     default-ON feature OFF;
//   * a negative value wrapped when cast to an unsigned type (=-1 became SIZE_MAX), and a
//     value past INT_MAX was undefined behaviour inside atoi;
//   * the first value read was kept for the whole process, so a driver that ran several
//     configurations in one process silently ran the first one every time.
//
// Every flag is now declared once, in the table in envflags.cpp, with its kind and its
// accepted range. main() validates the whole environment before any work starts
// (validateEnvironment) and prints one "[config] NAME=<parsed value>" line per flag that is
// set. The readers below parse the same way, per call, and never cache: a value changed
// between two calls in one process takes effect on the second call.
//
// Output contract: for every value the table accepts, a reader returns exactly what the
// release expression returned (atoi and strtoll agree on every base-10 integer inside int;
// atof is defined as strtod). A value it does not accept is a hard error -- exit status 2
// with a message naming the flag -- never a silent default. Only values that previously
// parsed to garbage, wrapped, or were undefined behaviour are rejected.
//
// Adding a flag: add one row to the table in envflags.cpp and read it with the reader that
// matches its kind, inside the function that uses it (never in a function-local static).
// tests/test_envflags.cpp scans src/ and fails if a TESSERACT_* literal is not in the
// table, if a reader does not match the table's kind, or if a static caches an env:: read.
// After the call-site retrofit, tests/test_envflags_retrofit.cpp also fails on any raw
// getenv of a flag and on any static that caches one.
#pragma once

#include <cstddef>
#include <cstdio>
#include <string>
#include <vector>

namespace ts { namespace env {

enum class Kind {
    Integer,   // base-10 integer, the whole string, inside [ilo, ihi]
    Real,      // decimal/exponent number (strtod), the whole string, finite, inside [rlo, rhi]
    Switch,    // base-10 integer inside int: 0 = off, any other value = on (release atoi != 0)
    Binary,    // exactly "0" (off) or "1" (on)
    Presence,  // set = on, whatever the value; diagnostic output only
    One,       // set = on; the only accepted value is "1" (set-means-on flags that change output)
    Text,      // free text: a path or a file name
    Choice,    // one word from a fixed list
    External   // read by the shipped wrapper scripts, never by this binary
};

struct Spec {
    const char* name;
    Kind kind;
    long long ilo, ihi;      // Integer bounds, inclusive
    double rlo, rhi;         // Real bounds, inclusive
    const char* choices;     // Choice: words separated by '|'
};

// The registry. `find` returns nullptr for a name that is not a TesserACT flag.
const Spec* find(const char* name);
const Spec* table(std::size_t& count);
const char* kindName(Kind k);

// A parsed value. `canonical` is what the [config] line prints.
struct Parsed {
    long long i = 0;
    double r = 0.0;
    bool on = false;
    std::string canonical;
};
// Pure parse of `value` against `spec`; no side effects. False with a reason on rejection.
bool parse(const Spec& spec, const char* value, Parsed& out, std::string& why);

// Readers. Each returns `def` when the flag is unset. A set value that the table rejects
// prints "error: NAME='VALUE' is not valid: ..." and exits with status 2. A name missing from
// the table, or read with the wrong reader for its kind, is a programming error: abort().
long long integer(const char* name, long long def);  // Kind::Integer
double real(const char* name, double def);           // Kind::Real
bool on(const char* name, bool def);                 // Kind::Switch or Kind::Binary
bool present(const char* name);                      // Kind::Presence or Kind::One
const char* text(const char* name);                  // Kind::Text or Kind::Choice; nullptr if unset
// build_v3: whether a registered flag is set at all (any kind). The fix switches use it to
// tell "unset, follow the umbrella" from an explicit value. It never parses the value; the
// matching reader does.
bool isSet(const char* name);

// Checks every TESSERACT_* variable in the environment. For each registered flag that is
// set, prints "[config] NAME=<parsed>" to `log`, sorted by name. A name that is not a
// TesserACT flag gets a warning line (it is ignored, as before). Returns 0 when every set
// flag is valid; otherwise prints one error line per invalid flag and returns 2.
int validateEnvironment(std::FILE* log);

// Phase 2 emit-B (F10 provenance): every TESSERACT_* variable of the environment, sorted by name,
// with its value, its parsed (canonical) value and kind when it is a registered flag. Reads the
// environment the way validateEnvironment does; no flag is interpreted here.
struct SetVariable {
    std::string name, value, canonical, kind;
    bool registered = false, valid = false;
};
std::vector<SetVariable> setVariables();

}}  // namespace ts::env
