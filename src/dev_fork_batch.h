// Developer-only, same-binary diagnostic batches. Never enabled by default.
#pragma once
#include <map>
#include <string>
#include <vector>
#include <cstdint>

namespace ts {
struct AssemblyOptions;
class UnitigGraph;
class SequenceStore;
namespace dev {
class ForkBatch {
public:
    enum class Entry { Disabled, Child, ParentComplete, Failed };
    bool prepare(const AssemblyOptions& options, std::string& error);
    Entry enter(AssemblyOptions& options, const UnitigGraph& graph,
                const SequenceStore& reads, double prefixSeconds, std::string& error);
    bool finishChild(double seconds, std::string& error);
    bool isChild() const { return child_; }
    ~ForkBatch();
private:
    struct Arm {
        std::string id;
        std::map<std::string, std::string> flags;
        long pid = 0;
        int exitCode = -1, signal = 0;
        uint64_t peakPss = 0;
        double seconds = 0, userSeconds = 0, systemSeconds = 0;
        std::string status = "pending";
    };
    bool writeBatch(const std::string& status, std::string& error) const;
    bool writeChild(const std::string& status, double seconds, std::string& error) const;
    bool enabled_ = false, child_ = false, childFinished_ = false;
    size_t childIndex_ = 0;
    long parentPid_ = 0;
    double prefixSeconds_ = 0;
    std::string manifest_, manifestHash_, root_, invocation_, binary_, binaryHash_, config_, state_, environment_;
    std::vector<Arm> arms_;
};
} // namespace dev
} // namespace ts
