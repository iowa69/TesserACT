// Phase 2, emit-B: a minimal JSON object builder for the report.json "p2" block.
#pragma once

#include <cstdint>
#include <cstdio>
#include <cmath>
#include <string>
#include <vector>

namespace ts {
namespace p2 {

inline std::string jsonEscape(const std::string& s) {
    std::string o;
    o.reserve(s.size() + 2);
    for (unsigned char c : s) {
        switch (c) {
            case '"': o += "\\\""; break;
            case '\\': o += "\\\\"; break;
            case '\n': o += "\\n"; break;
            case '\r': o += "\\r"; break;
            case '\t': o += "\\t"; break;
            default:
                if (c < 0x20) {
                    char b[8];
                    std::snprintf(b, sizeof b, "\\u%04x", c);
                    o += b;
                } else {
                    o += static_cast<char>(c);
                }
        }
    }
    return o;
}

// {"a": 1, "b": "x", "c": {...}} built field by field; values are written compactly.
class JsonObj {
public:
    JsonObj& raw(const std::string& key, const std::string& value) {
        if (!body_.empty()) body_ += ", ";
        body_ += "\"" + jsonEscape(key) + "\": " + value;
        return *this;
    }
    JsonObj& str(const std::string& key, const std::string& v) { return raw(key, "\"" + jsonEscape(v) + "\""); }
    JsonObj& u(const std::string& key, uint64_t v) { return raw(key, std::to_string(v)); }
    JsonObj& i(const std::string& key, long long v) { return raw(key, std::to_string(v)); }
    JsonObj& b(const std::string& key, bool v) { return raw(key, v ? "true" : "false"); }
    JsonObj& num(const std::string& key, double v, int prec = 6) {
        if (!std::isfinite(v)) return raw(key, "null");
        char buf[64];
        std::snprintf(buf, sizeof buf, "%.*f", prec, v);
        return raw(key, buf);
    }
    JsonObj& strs(const std::string& key, const std::vector<std::string>& v) {
        std::string a = "[";
        for (size_t j = 0; j < v.size(); ++j) a += (j ? ", \"" : "\"") + jsonEscape(v[j]) + "\"";
        return raw(key, a + "]");
    }
    std::string done() const { return "{" + body_ + "}"; }
    bool empty() const { return body_.empty(); }

private:
    std::string body_;
};

inline std::string jsonArray(const std::vector<std::string>& items) {
    std::string a = "[";
    for (size_t j = 0; j < items.size(); ++j) a += (j ? ", " : "") + items[j];
    return a + "]";
}

}  // namespace p2
}  // namespace ts
