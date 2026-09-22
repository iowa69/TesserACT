// Deliberately include the implementation to exercise its private hash utility;
// link this fixture without dev_fork_batch.o. No test hooks enter production.
#include "../src/dev_fork_batch.cpp"
#include <atomic>
#include <thread>

int main(int argc, char** argv) {
    if (argc == 3 && std::string(argv[1]) == "hash") {
        std::string hash;
        if (!ts::dev::fileHash(argv[2], hash)) return 2;
        std::puts(hash.c_str()); return 0;
    }
    if (argc != 5) return 2;
    const std::string action = argv[1];
    ts::AssemblyOptions options;
    options.outDir = argv[2]; options.commandLine = "developer fixture " + action;
    options.threads = 1;
    options.maxMemoryBytes = 9007199254740993ULL; // Exact integer provenance beyond double precision.
    options.libraries.push_back({argv[3], argv[4]});
    ts::SequenceStore reads; std::string error;
    if (!reads.load(options.libraries, 1, error)) return 2;
    ts::UnitigGraph graph; graph.setK(31);
    graph.nodes.push_back({std::string(80, 'A'), 10, false, {}});
    ts::dev::ForkBatch batch;
    if (!batch.prepare(options, error)) { std::fprintf(stderr, "%s\n", error.c_str()); return 3; }
    if (action == "prepare_hold") {
        std::ofstream(options.outDir + "/claim_ready") << "claimed\n";
        for (;;) ::pause();
    }
    std::atomic<bool> stop{false}; std::thread other;
    if (action == "thread") other = std::thread([&] { while (!stop) std::this_thread::yield(); });
    if (action == "sigchld_ignored") ::signal(SIGCHLD, SIG_IGN);
    if (action == "sigchld_no_wait") {
        struct sigaction disposition{}; disposition.sa_handler = SIG_DFL;
        disposition.sa_flags = SA_NOCLDWAIT; sigemptyset(&disposition.sa_mask);
        ::sigaction(SIGCHLD, &disposition, nullptr);
    }
    if (action == "blocked") {
        sigset_t mask; sigemptyset(&mask); sigaddset(&mask, SIGINT); sigaddset(&mask, SIGTERM);
        ::sigprocmask(SIG_BLOCK, &mask, nullptr);
    }
    if (action == "ignored_cancel") { ::signal(SIGINT, SIG_IGN); ::signal(SIGTERM, SIG_IGN); }
    auto entry = batch.enter(options, graph, reads, 12.5, error);
    if (other.joinable()) { stop = true; other.join(); }
    if (entry == ts::dev::ForkBatch::Entry::Failed) { std::fprintf(stderr, "%s\n", error.c_str()); return 1; }
    if (entry == ts::dev::ForkBatch::Entry::ParentComplete) return 0;
    if (entry != ts::dev::ForkBatch::Entry::Child) return 4;
    if (action == "exit") return 7;
    if (action == "signal") { ::raise(SIGTERM); return 9; }
    if (action == "cancel" || action == "ignored_cancel") { for (;;) ::pause(); }
    if (action != "success") return 5;
    // Dirty private pages long enough for the parent's memory sampler to see.
    std::vector<unsigned char> scratch(16 * 1024 * 1024, 137);
    for (size_t i = 0; i < scratch.size(); i += 4096) scratch[i] = 19;
    struct timespec delay{0, 120000000}; ::nanosleep(&delay, nullptr);
    std::ofstream(options.outDir + "/contigs.fasta") << ">synthetic\nACGT\n";
    std::ofstream report(options.outDir + "/report.json");
    report << "{";
    bool first = true;
    for (const char* flag : ts::dev::flags) {
        if (!first) report << ',';
        first = false;
        const char* value = std::getenv(flag);
        report << ts::dev::quote(flag) << ':' << (value ? ts::dev::quote(value) : "null");
    }
    report << "}\n"; report.close();
    return batch.finishChild(.125, error) ? 0 : 6;
}
