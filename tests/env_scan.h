// Test helper (header only; not a test of its own): a source scan for TESSERACT_* flag use.
//
// Always checked:
//   * every "TESSERACT_*" string literal in src/ is in the flag table (src/envflags.cpp);
//   * every env::integer/real/on/present/text("NAME") call (and replicon's envSize) uses the
//     reader that matches NAME's kind in the table;
//   * no static caches a value read through the env:: readers;
//   * the flags dev_fork_batch varies between children are registered and read "0"/"1".
// Strict (after the phase-3 retrofit):
//   * no raw getenv outside envflags.cpp (except HOME in main.cpp);
//   * no static caches a getenv result either;
//   * every registered flag (other than wrapper-script variables) is read somewhere.
#pragma once

#include <algorithm>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <map>
#include <regex>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "envflags.h"

namespace envscan {

inline std::string readFile(const std::filesystem::path& p) {
    std::ifstream in(p, std::ios::binary);
    std::ostringstream s;
    s << in.rdbuf();
    return s.str();
}

// Blank out comments; with `strings`, blank string and char literal contents too. Offsets
// and newlines are kept so positions still map to lines.
inline std::string blank(const std::string& src, bool strings) {
    std::string out = src;
    enum { Code, Line, Block, Str, Chr } st = Code;
    for (size_t i = 0; i < src.size(); ++i) {
        const char c = src[i], d = i + 1 < src.size() ? src[i + 1] : '\0';
        switch (st) {
            case Code:
                if (c == '/' && d == '/') { st = Line; out[i] = ' '; }
                else if (c == '/' && d == '*') { st = Block; out[i] = ' '; }
                else if (c == '"') st = Str;
                else if (c == '\'') st = Chr;
                break;
            case Line:
                if (c == '\n') st = Code; else out[i] = ' ';
                break;
            case Block:
                if (c == '*' && d == '/') { out[i] = ' '; out[i + 1] = ' '; ++i; st = Code; }
                else if (c != '\n') out[i] = ' ';
                break;
            case Str:
            case Chr:
                if (c == '\\') { if (strings) { out[i] = ' '; if (i + 1 < src.size()) out[i + 1] = ' '; } ++i; }
                else if ((st == Str && c == '"') || (st == Chr && c == '\'')) st = Code;
                else if (strings && c != '\n') out[i] = ' ';
                break;
        }
    }
    return out;
}

inline int lineOf(const std::string& s, size_t pos) {
    return 1 + static_cast<int>(std::count(s.begin(), s.begin() + static_cast<long>(pos), '\n'));
}

inline bool isIdent(char c) { return std::isalnum(static_cast<unsigned char>(c)) || c == '_'; }

// Violations found in the sources directly under `dir`, one line each.
inline std::vector<std::string> scan(const std::filesystem::path& dir, bool strict) {
    std::vector<std::string> v;
    std::vector<std::filesystem::path> files;
    if (std::filesystem::is_directory(dir)) {
        for (const auto& e : std::filesystem::directory_iterator(dir)) {
            const std::string ext = e.path().extension().string();
            if (ext == ".cpp" || ext == ".h") files.push_back(e.path());
        }
    }
    std::sort(files.begin(), files.end());
    if (files.empty()) { v.push_back("no sources under " + dir.string()); return v; }

    const std::regex literal("\"(TESSERACT_[A-Z0-9_]+)\"");
    const std::regex reader("\\b(env::integer|env::real|env::on|env::present|env::text|envSize)\\s*\\(\\s*\"(TESSERACT_[A-Z0-9_]+)\"");
    const std::regex rawGetenv("getenv\\s*\\(");
    const std::regex indirect(
        "\\b(fixflags::fixLevel|fixEnabled)\\s*\\(\\s*\"(TESSERACT_[A-Z0-9_]+)\"|"
        "constexpr const char\\* k[A-Za-z0-9]+\\s*=\\s*\"(TESSERACT_[A-Z0-9_]+)\"");
    using ts::env::Kind;
    const std::map<std::string, std::set<Kind>> readerKinds = {
        {"env::integer", {Kind::Integer}}, {"envSize", {Kind::Integer}},
        {"env::real", {Kind::Real}},
        {"env::on", {Kind::Switch, Kind::Binary}},
        {"env::present", {Kind::Presence, Kind::One}},
        {"env::text", {Kind::Text, Kind::Choice}}};
    std::set<std::string> readNames;

    for (const auto& f : files) {
        const std::string name = f.filename().string();
        const std::string raw = readFile(f);
        const std::string code = blank(raw, false);   // comments gone, literals kept
        const std::string bare = blank(raw, true);    // comments and literal contents gone
        const bool isRegistry = name == "envflags.cpp" || name == "envflags.h";
        const std::string at = name + ":";

        for (auto it = std::sregex_iterator(code.begin(), code.end(), literal); it != std::sregex_iterator(); ++it) {
            const std::string flag = (*it)[1];
            if (!ts::env::find(flag.c_str()))
                v.push_back(at + std::to_string(lineOf(code, static_cast<size_t>(it->position()))) + ": " + flag +
                            " is not in the flag table (src/envflags.cpp)");
        }
        for (auto it = std::sregex_iterator(code.begin(), code.end(), reader); it != std::sregex_iterator(); ++it) {
            const std::string fn = (*it)[1], flag = (*it)[2];
            readNames.insert(flag);
            const ts::env::Spec* s = ts::env::find(flag.c_str());
            if (s && !readerKinds.at(fn).count(s->kind))
                v.push_back(at + std::to_string(lineOf(code, static_cast<size_t>(it->position()))) + ": " + flag +
                            " is a " + ts::env::kindName(s->kind) + " flag read through " + fn);
        }
        // build_v3: the defect-fix switches are read through umbrella-aware helpers that call
        // env::isSet/on/integer with the name they are given -- fixflags::fixLevel("NAME")
        // (graph_fix_flags.h), fixEnabled("NAME") (resolve.cpp), and the name constants of
        // emit_fixflags.h. Those count as reads, and the flag must be Binary or Integer.
        for (auto it = std::sregex_iterator(code.begin(), code.end(), indirect); it != std::sregex_iterator(); ++it) {
            const std::string flag = (*it)[2].matched ? std::string((*it)[2]) : std::string((*it)[3]);
            readNames.insert(flag);
            const ts::env::Spec* s = ts::env::find(flag.c_str());
            if (s && s->kind != Kind::Binary && s->kind != Kind::Integer)
                v.push_back(at + std::to_string(lineOf(code, static_cast<size_t>(it->position()))) + ": " + flag +
                            " is a " + ts::env::kindName(s->kind) + " flag read through a fix-switch helper");
        }
        if (isRegistry) continue;
        if (strict) {
            // Raw getenv: only envflags.cpp may call it, plus HOME in main.cpp.
            for (auto it = std::sregex_iterator(code.begin(), code.end(), rawGetenv); it != std::sregex_iterator(); ++it) {
                const size_t pos = static_cast<size_t>(it->position() + it->length());
                if (code.compare(pos, 6, "\"HOME\"") == 0) continue;
                v.push_back(at + std::to_string(lineOf(code, static_cast<size_t>(it->position()))) +
                            ": raw getenv (use the env:: readers in envflags.h)");
            }
        }
        // A static whose initialiser reads the environment caches the first value for the
        // whole process. Find each `static` declaration's initialiser and look inside it.
        for (size_t p = bare.find("static"); p != std::string::npos; p = bare.find("static", p + 6)) {
            if ((p > 0 && isIdent(bare[p - 1])) || (p + 6 < bare.size() && isIdent(bare[p + 6]))) continue;
            int depth = 0;
            bool init = false;
            size_t end = p + 6, start = std::string::npos;
            for (; end < bare.size(); ++end) {
                const char c = bare[end];
                if (c == '(' || c == '[' || c == '{') {
                    if (c == '{' && depth == 0 && !init) {
                        // `static T f(...) {` or `static struct X {` is a body, not a value;
                        // `static T x{...}` is brace initialisation.
                        size_t q = end;
                        while (q > p && std::isspace(static_cast<unsigned char>(bare[q - 1]))) --q;
                        const std::string before = bare.substr(p, q - p);
                        const bool body = bare[q - 1] == ')' || before.find('(') != std::string::npos ||
                                          before.find("struct") != std::string::npos ||
                                          before.find("class") != std::string::npos ||
                                          before.find("namespace") != std::string::npos;
                        if (body) break;
                        init = true;
                        start = end;
                    }
                    ++depth;
                } else if (c == ')' || c == ']' || c == '}') {
                    --depth;
                } else if (c == '=' && depth == 0 && !init) {
                    init = true;
                    start = end;
                } else if (c == ';' && depth == 0) {
                    break;
                }
            }
            if (!init || start == std::string::npos) continue;
            const std::string initialiser = bare.substr(start, end - start);
            const bool readsEnv = initialiser.find("env::") != std::string::npos ||
                                  (strict && initialiser.find("getenv") != std::string::npos);
            if (readsEnv)
                v.push_back(at + std::to_string(lineOf(bare, p)) +
                            ": static caches an environment flag (read it per call instead)");
        }
    }

    size_t n = 0;
    const ts::env::Spec* table = ts::env::table(n);
    if (strict) {
        for (size_t i = 0; i < n; ++i) {
            if (table[i].kind == Kind::External) continue;
            if (!readNames.count(table[i].name))
                v.push_back(std::string("envflags.cpp: ") + table[i].name + " is registered but never read");
        }
    }
    // The developer fork batch varies flags between children after fork(). They must be
    // registered flags that read "0"/"1" as off/on (the batch writes nothing else).
    if (std::filesystem::exists(dir / "dev_fork_batch.cpp")) {
        const std::string fb = blank(readFile(dir / "dev_fork_batch.cpp"), false);
        const size_t a = fb.find("flags = {{");
        const size_t b = a == std::string::npos ? std::string::npos : fb.find("}};", a);
        if (a == std::string::npos || b == std::string::npos) v.push_back("dev_fork_batch.cpp: flag list not found");
        else {
            const std::string list = fb.substr(a, b - a);
            int count = 0;
            for (auto it = std::sregex_iterator(list.begin(), list.end(), literal); it != std::sregex_iterator(); ++it) {
                ++count;
                const ts::env::Spec* s = ts::env::find(std::string((*it)[1]).c_str());
                if (!s || (s->kind != Kind::Binary && s->kind != Kind::Switch))
                    v.push_back("dev_fork_batch.cpp: " + std::string((*it)[1]) + " must be a registered 0/1 flag");
            }
            if (count == 0) v.push_back("dev_fork_batch.cpp: empty flag list");
        }
    }
    return v;
}

// The sources this test was built from: tests/ sits next to src/, and make runs the
// component tests from the tree root.
inline std::filesystem::path sourceDir(const char* thisFile) {
    return std::filesystem::path(thisFile).parent_path().parent_path() / "src";
}

}  // namespace envscan
