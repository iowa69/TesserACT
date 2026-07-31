#include "util.h"

#include <algorithm>
#include <cstdio>
#include <sys/resource.h>
#include <sys/stat.h>

namespace ts {
namespace util {

std::string commify(long long v) {
    bool neg = v < 0;
    unsigned long long x = neg ? static_cast<unsigned long long>(-(v + 1)) + 1
                               : static_cast<unsigned long long>(v);
    std::string digits = std::to_string(x);
    std::string out;
    int count = 0;
    for (int i = static_cast<int>(digits.size()) - 1; i >= 0; --i) {
        out += digits[static_cast<size_t>(i)];
        if (++count % 3 == 0 && i > 0) out += ',';
    }
    if (neg) out += '-';
    std::reverse(out.begin(), out.end());
    return out;
}

std::string humanBytes(double bytes) {
    static const char* units[] = {"B", "KB", "MB", "GB", "TB"};
    int u = 0;
    while (bytes >= 1024.0 && u < 4) { bytes /= 1024.0; ++u; }
    char buf[64];
    std::snprintf(buf, sizeof(buf), u == 0 ? "%.0f %s" : "%.2f %s", bytes, units[u]);
    return buf;
}

bool fileExists(const std::string& path) {
    struct stat st;
    return ::stat(path.c_str(), &st) == 0;
}

bool makeDirs(const std::string& path) {
    if (path.empty() || path == ".") return true;
    if (fileExists(path)) return true;
    std::string acc;
    for (size_t i = 0; i < path.size(); ++i) {
        acc += path[i];
        if (path[i] == '/' || i + 1 == path.size()) {
            if (acc == "/" || acc.empty()) continue;
            if (!fileExists(acc) && ::mkdir(acc.c_str(), 0755) != 0) return false;
        }
    }
    return fileExists(path);
}

std::vector<std::string> split(const std::string& s, char sep) {
    std::vector<std::string> out;
    size_t start = 0;
    while (true) {
        size_t pos = s.find(sep, start);
        if (pos == std::string::npos) { out.push_back(s.substr(start)); break; }
        out.push_back(s.substr(start, pos - start));
        start = pos + 1;
    }
    return out;
}

std::string dirname(const std::string& path) {
    size_t pos = path.find_last_of('/');
    return pos == std::string::npos ? "." : path.substr(0, pos);
}

long long peakMemoryBytes() {
    struct rusage ru;
    if (getrusage(RUSAGE_SELF, &ru) != 0) return 0;
    // ru_maxrss is in kilobytes on Linux.
    return static_cast<long long>(ru.ru_maxrss) * 1024;
}

}  // namespace util
}  // namespace ts
