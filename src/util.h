// Small shared helpers.
#pragma once

#include <chrono>
#include <string>
#include <vector>

namespace ts {
namespace util {

class Timer {
public:
    Timer() : start_(std::chrono::steady_clock::now()) {}
    void reset() { start_ = std::chrono::steady_clock::now(); }
    double elapsed() const {
        return std::chrono::duration<double>(std::chrono::steady_clock::now() - start_).count();
    }
private:
    std::chrono::steady_clock::time_point start_;
};

std::string commify(long long v);
std::string humanBytes(double bytes);
bool makeDirs(const std::string& path);
bool fileExists(const std::string& path);
std::vector<std::string> split(const std::string& s, char sep);
std::string dirname(const std::string& path);

// Peak resident set size in bytes, or 0 when unavailable.
long long peakMemoryBytes();

}  // namespace util
}  // namespace ts
