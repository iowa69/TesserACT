#include "dev_fork_batch.h"
#include "assembler.h"
#include "util.h"

#include <algorithm>
#include <array>
#include <cerrno>
#include <chrono>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <dirent.h>
#include <fcntl.h>
#include <fstream>
#include <iomanip>
#include <limits.h>
#include <sstream>
#include <sys/prctl.h>
#include <sys/resource.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <signal.h>
#include <time.h>
#include <unistd.h>

namespace ts { namespace dev { namespace {
const std::array<const char*, 7> flags = {{
    "TESSERACT_ROUTE_DISTANCE", "TESSERACT_WEIGHTED_RESOLVER_COVERAGE",
    "TESSERACT_WEIGHTED_ELIGIBLE_COVERAGE", "TESSERACT_EXACT_READ_THREADS",
    "TESSERACT_OWNED_ANCHOR_PREFIX", "TESSERACT_EXCLUDE_SHARED_REPEAT_SUPPORT",
    "TESSERACT_OBSERVED_GRAPH_COVERAGE"
}};

std::string quote(const std::string& value) {
    std::string out = "\"";
    for (unsigned char c : value) {
        if (c == '\\' || c == '"') { out += '\\'; out += char(c); }
        else if (c < 32) { char b[7]; std::snprintf(b, sizeof b, "\\u%04x", unsigned(c)); out += b; }
        else out += char(c);
    }
    return out + '"';
}
template<class T> std::string number(T value) {
    std::ostringstream out; out.imbue(std::locale::classic());
    out << std::setprecision(17) << value; return out.str();
}
std::string absolute(const std::string& path) {
    char buf[PATH_MAX];
    return ::realpath(path.c_str(), buf) ? std::string(buf) : std::string();
}
bool exists(const std::string& path) { struct stat st; return ::lstat(path.c_str(), &st) == 0; }

// Small streaming SHA-256 implementation for executable/manifest identities.
// This is provenance hashing, not a substitute for the caller's source manifest.
class Sha256 {
    uint32_t h[8] = {0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19};
    uint64_t bytes = 0; size_t used = 0; unsigned char block[64]{};
    static uint32_t rotr(uint32_t x, int n) { return (x >> n) | (x << (32-n)); }
    void compress() {
        static const uint32_t k[64] = {
            0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
            0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
            0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
            0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
            0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
            0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
            0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
            0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2};
        uint32_t w[64];
        for (int i=0;i<16;++i) w[i]=(uint32_t(block[4*i])<<24)|(uint32_t(block[4*i+1])<<16)|(uint32_t(block[4*i+2])<<8)|block[4*i+3];
        for (int i=16;i<64;++i) {
            uint32_t a=w[i-15],b=w[i-2];
            w[i]=w[i-16]+(rotr(a,7)^rotr(a,18)^(a>>3))+w[i-7]+(rotr(b,17)^rotr(b,19)^(b>>10));
        }
        uint32_t a=h[0],b=h[1],c=h[2],d=h[3],e=h[4],f=h[5],g=h[6],z=h[7];
        for (int i=0;i<64;++i) {
            uint32_t t1=z+(rotr(e,6)^rotr(e,11)^rotr(e,25))+((e&f)^(~e&g))+k[i]+w[i];
            uint32_t t2=(rotr(a,2)^rotr(a,13)^rotr(a,22))+((a&b)^(a&c)^(b&c));
            z=g;g=f;f=e;e=d+t1;d=c;c=b;b=a;a=t1+t2;
        }
        h[0]+=a;h[1]+=b;h[2]+=c;h[3]+=d;h[4]+=e;h[5]+=f;h[6]+=g;h[7]+=z;
    }
public:
    void add(const void* p, size_t n) {
        const auto* s=static_cast<const unsigned char*>(p); bytes+=n;
        while(n) { size_t take=std::min(n,64-used); std::memcpy(block+used,s,take);used+=take;s+=take;n-=take;if(used==64){compress();used=0;} }
    }
    std::string finish() {
        uint64_t bits=bytes*8; unsigned char one=0x80,zero=0;add(&one,1);
        while(used!=56)add(&zero,1);
        unsigned char end[8];for(int i=0;i<8;++i)end[7-i]=static_cast<unsigned char>(bits>>(i*8));add(end,8);
        char out[65];for(int i=0;i<8;++i)std::snprintf(out+i*8,9,"%08x",h[i]);return out;
    }
};
bool fileHash(const std::string& path, std::string& result) {
    std::ifstream input(path,std::ios::binary); if(!input)return false;
    Sha256 hash;char buffer[65536];
    while(input){input.read(buffer,sizeof buffer);hash.add(buffer,static_cast<size_t>(input.gcount()));}
    if (!input.eof()) return false;
    result = hash.finish();
    return true;
}
bool atomicWrite(const std::string& path,const std::string& text,std::string& error) {
    const std::string tmp=path+".tmp."+std::to_string(::getpid());
    std::FILE* f=std::fopen(tmp.c_str(),"wb");
    if(!f){error="cannot write diagnostic provenance: "+tmp;return false;}
    bool ok=std::fwrite(text.data(),1,text.size(),f)==text.size();
    if(std::fclose(f)!=0)ok=false;
    if(!ok || ::rename(tmp.c_str(),path.c_str())!=0){::unlink(tmp.c_str());error="cannot complete diagnostic provenance: "+path;return false;}
    return true;
}
std::string optionJson(const AssemblyOptions& o) {
    std::map<std::string,std::string> m;
#define STR(x) m[#x]=quote(o.x)
#define NUM(x) m[#x]=number(o.x)
#define BOOL(x) m[#x]=o.x?"true":"false"
    STR(outDir);STR(commandLine);STR(organism);STR(organismModelPath);STR(isPanelPath);STR(isSitesPath);STR(qcPath);STR(mapperDir);
    NUM(threads);NUM(forcedCutoff);NUM(trustCutoff);NUM(minMaskRun);NUM(minContigLen);NUM(minLinkSupport);NUM(linkSupportPerX);NUM(tieRatio);
    NUM(bubbleCoverageLimit);NUM(simplifyRounds);NUM(polishPasses);NUM(maxMemoryBytes);NUM(ladderUnionMaxPresent);NUM(ladderUnionMinLen);NUM(trimTerminalOverlap);
    BOOL(verbose);BOOL(correctReads);BOOL(resolveRepeats);BOOL(scaffold);BOOL(gapFill);BOOL(polish);BOOL(emitGfa);BOOL(emitHtml);BOOL(emitUnitigs);
    BOOL(ladderUnion);BOOL(dedupContained);BOOL(layout);BOOL(userSetK);BOOL(userSetMinLink);BOOL(userSetTie);BOOL(userSetLinkPerX);BOOL(userSetBubble);BOOL(userSetRounds);BOOL(userSetPolishPasses);
#undef STR
#undef NUM
#undef BOOL
    m["mode"]=quote(runModeName(o.mode));m["mapPolisher"]=quote(mapperName(o.mapPolisher));
    m["qtrim"]="{\"enabled\":"+std::string(o.qtrim.enabled?"true":"false")+",\"window\":"+std::to_string(o.qtrim.windowSize)+",\"quality\":"+std::to_string(o.qtrim.meanQuality)+",\"offset\":"+std::to_string(o.qtrim.phredOffset)+"}";
    std::string ks="[";for(int k:o.kValues){if(ks.size()>1)ks+=',';ks+=std::to_string(k);}m["kValues"]=ks+"]";
    std::string libs="[";for(const auto& l:o.libraries){if(libs.size()>1)libs+=',';libs+="{\"r1\":"+quote(l.r1)+",\"r2\":"+quote(l.r2)+",\"interleaved\":"+(l.interleaved?"true":"false")+",\"insertMean\":"+number(l.insertMean)+",\"insertStdDev\":"+number(l.insertStdDev)+",\"readLength\":"+std::to_string(l.readLength)+",\"oriented\":"+(l.oriented?"true":"false")+"}";}m["libraries"]=libs+"]";
    std::string out="{";for(const auto& p:m){if(out.size()>1)out+=',';out+=quote(p.first)+":"+p.second;}return out+"}";
}
std::string flagJson(const std::map<std::string,std::string>& values) {
    std::string out="{";for(const char* f:flags){if(out.size()>1)out+=',';auto p=values.find(f);out+=quote(f)+":"+(p==values.end()?"null":quote(p->second));}return out+"}";
}
size_t taskCount() {
    DIR* d=::opendir("/proc/self/task");if(!d)return 0;size_t n=0;
    while(auto* e=::readdir(d)){if(e->d_name[0]>='0'&&e->d_name[0]<='9')++n;}::closedir(d);return n;
}
uint64_t pss(long pid) {
    std::ifstream f("/proc/"+std::to_string(pid)+"/smaps_rollup");std::string line;
    while(std::getline(f,line)){if(line.rfind("Pss:",0)==0){std::istringstream in(line.substr(4));uint64_t kb=0;in>>kb;return kb*1024;}}return 0;
}
volatile sig_atomic_t interrupted=0;
// Only record here. Forwarding after wait4 has reaped a child could target a
// reused PID; the synchronous wait loop forwards only while still unreaped.
void interruptedHandler(int signal) { interrupted=signal; }
bool safeSignalState(std::string& error) {
    struct sigaction childAction{};
    if (::sigaction(SIGCHLD, nullptr, &childAction) != 0 ||
        childAction.sa_handler != SIG_DFL || (childAction.sa_flags & SA_NOCLDWAIT)) {
        error = "developer fork requires default, waitable SIGCHLD handling";
        return false;
    }
    sigset_t mask;
    if (::sigprocmask(SIG_SETMASK, nullptr, &mask) != 0 ||
        ::sigismember(&mask, SIGINT) || ::sigismember(&mask, SIGTERM)) {
        error = "developer fork requires unblocked SIGINT and SIGTERM";
        return false;
    }
    return true;
}
class SignalScope {
    struct sigaction oldInt{},oldTerm{};bool ready=false;
public:
    bool install(){struct sigaction action{};action.sa_handler=interruptedHandler;sigemptyset(&action.sa_mask);if(sigaction(SIGINT,&action,&oldInt)!=0)return false;if(sigaction(SIGTERM,&action,&oldTerm)!=0){sigaction(SIGINT,&oldInt,nullptr);return false;}ready=true;return true;}
    void restore(){if(ready){sigaction(SIGINT,&oldInt,nullptr);sigaction(SIGTERM,&oldTerm,nullptr);ready=false;}}
    ~SignalScope(){restore();}
};
} // namespace

bool ForkBatch::prepare(const AssemblyOptions& o,std::string& error) {
    const char* requested=std::getenv("TESSERACT_DEV_FORK_BATCH");if(!requested)return true;
    enabled_=true;
    if(!*requested){error="TESSERACT_DEV_FORK_BATCH requires a manifest path";return false;}
    if(!o.organismModelPath.empty()||!o.organism.empty()||!o.isPanelPath.empty()||!o.isSitesPath.empty()||o.mapPolisher!=Mapper::None||!o.resolveRepeats){error="developer fork batch requires model-free paired resolution and no external mapper";return false;}
    if(std::getenv("TESSERACT_PLASMID_CLUSTERS")||std::getenv("TESSERACT_JOIN_DUMP")){error="developer fork batch forbids external dump/model-cluster paths";return false;}
    manifest_=absolute(requested);root_=absolute(o.outDir);invocation_=o.commandLine;config_=optionJson(o);
    if(manifest_.empty()||root_.empty()){error="cannot resolve developer manifest/output directory";return false;}
    std::ifstream manifestFile(manifest_, std::ios::binary);
    std::string manifestText; char chunk[4096];
    while (manifestFile) {
        manifestFile.read(chunk, sizeof chunk);
        manifestText.append(chunk, static_cast<size_t>(manifestFile.gcount()));
        if (manifestText.size() > 1048576) { error = "developer manifest exceeds 1 MiB"; return false; }
    }
    if (!manifestFile.eof()) { error = "cannot read developer manifest"; return false; }
    Sha256 manifestDigest; manifestDigest.add(manifestText.data(), manifestText.size());
    manifestHash_ = manifestDigest.finish();
    std::istringstream input(manifestText);std::string line;if(!std::getline(input,line)||line!="arm\tflags"){error="developer manifest must start with arm<TAB>flags";return false;}
    size_t lines=1,bytes=line.size();
    while(std::getline(input,line)){
        ++lines;bytes+=line.size();if(bytes>1048576||lines>129){error="developer manifest exceeds 128 arms or 1 MiB";return false;}
        size_t tab=line.find('\t');if(tab==std::string::npos||line.find('\t',tab+1)!=std::string::npos){error="developer manifest needs exactly two tab-separated columns";return false;}
        Arm a;a.id=line.substr(0,tab);std::string values=line.substr(tab+1);
        if(a.id.empty()||a.id.size()>80||!std::isalnum(static_cast<unsigned char>(a.id[0]))||!std::all_of(a.id.begin(),a.id.end(),[](unsigned char c){return std::isalnum(c)||c=='_'||c=='-';})){error="invalid developer arm identifier";return false;}
        if(std::any_of(arms_.begin(),arms_.end(),[&](const Arm& b){return a.id==b.id;})){error="duplicate developer arm identifier";return false;}
        if(values!="-"){
            size_t start=0;
            while(true){size_t comma=values.find(',',start);std::string token=values.substr(start,comma==std::string::npos?comma:comma-start);size_t eq=token.find('=');
                std::string key=token.substr(0,eq),value=eq==std::string::npos?"":token.substr(eq+1);
                if(std::find(flags.begin(),flags.end(),key)==flags.end()||(value!="0"&&value!="1")||!a.flags.emplace(key,value).second){error="forbidden, duplicate or invalid developer flag: "+token;return false;}
                if (comma == std::string::npos) break;
                start = comma + 1;
            }
        }
        if(exists(root_+"/"+a.id)){error="developer arm output already exists: "+a.id;return false;}
        arms_.push_back(std::move(a));
    }
    if(!input.eof()||arms_.empty()){error="developer manifest is empty or unreadable";return false;}
    for(const char* old:{"contigs.fasta","scaffolds.fasta","assembly_graph.gfa","report.json","fork_batch.json"})if(exists(root_+"/"+old)){error="developer batch output root contains existing assembly/batch output";return false;}
    char path[PATH_MAX];ssize_t n=::readlink("/proc/self/exe",path,sizeof(path)-1);
    if(n<=0){error="cannot identify developer batch executable";return false;}path[n]=0;binary_=path;
    if(!fileHash("/proc/self/exe",binaryHash_)){error="cannot hash developer executable/manifest";return false;}
    std::map<std::string, std::string> env;
    for (char** entry = ::environ; entry && *entry; ++entry) {
        std::string pair(*entry);
        if (pair.rfind("TESSERACT_", 0) != 0) continue;
        const size_t eq = pair.find('=');
        if (eq != std::string::npos) env[pair.substr(0, eq)] = pair.substr(eq+1);
    }
    environment_ = "{";
    for (const auto& e : env) {
        if (environment_.size() > 1) environment_ += ',';
        environment_ += quote(e.first) + ":" + quote(e.second);
    }
    environment_ += "}";
    parentPid_ = ::getpid();
    // Persistent exclusive claim is acquired before preprocessing. Never remove
    // it on failure: a partial root must be explicitly replaced by a fresh one.
    const std::string claim = root_ + "/.fork_batch_claim";
    int claimFd = ::open(claim.c_str(), O_CREAT | O_EXCL | O_WRONLY | O_CLOEXEC, 0644);
    if (claimFd < 0) { error = "developer batch output root is already claimed or unwritable"; return false; }
    const std::string owner = "{\"parent_pid\":" + std::to_string(parentPid_) +
        ",\"actual_parent_invocation\":" + quote(invocation_) + "}\n";
    const bool written = ::write(claimFd, owner.data(), owner.size()) == static_cast<ssize_t>(owner.size());
    const bool closed = ::close(claimFd) == 0;
    if (!written || !closed) { error = "cannot record developer batch output claim"; return false; }
    return true;
}

bool ForkBatch::writeBatch(const std::string& status,std::string& error) const {
    std::string text="{\n\"schema\":1,\"mode\":\"fork-batch\",\"diagnostic_only\":true,\"status\":"+quote(status)+",\"parent_pid\":"+std::to_string(parentPid_)+",\"actual_parent_invocation\":"+quote(invocation_)+",\"manifest\":"+quote(manifest_)+",\"manifest_sha256\":"+quote(manifestHash_)+",\"binary\":"+quote(binary_)+",\"binary_sha256\":"+quote(binaryHash_)+",\"shared_prefix_seconds\":"+number(prefixSeconds_)+",\"effective_parent_options\":"+config_+",\"parent_tesseract_environment\":"+environment_+",\"shared_state\":"+state_+",\"arms\":[";
    for(size_t i=0;i<arms_.size();++i){const auto& a=arms_[i];if(i)text+=',';text+="{\"arm\":"+quote(a.id)+",\"flags\":"+flagJson(a.flags)+",\"pid\":"+std::to_string(a.pid)+",\"status\":"+quote(a.status)+",\"exit_code\":"+std::to_string(a.exitCode)+",\"signal\":"+std::to_string(a.signal)+",\"waited_wall_seconds\":"+number(a.seconds)+",\"child_user_seconds\":"+number(a.userSeconds)+",\"child_system_seconds\":"+number(a.systemSeconds)+",\"sampled_parent_child_pss_sum_peak_bytes\":"+std::to_string(a.peakPss)+"}";}
    text+="]}\n";return atomicWrite(root_+"/fork_batch.json",text,error);
}
bool ForkBatch::writeChild(const std::string& status,double seconds,std::string& error) const {
    const auto& a=arms_[childIndex_];
    std::string text="{\n\"schema\":1,\"mode\":\"fork-child\",\"diagnostic_only\":true,\"status\":"+quote(status)+",\"arm\":"+quote(a.id)+",\"parent_pid\":"+std::to_string(parentPid_)+",\"pid\":"+std::to_string(::getpid())+",\"batch_provenance\":"+quote(root_+"/fork_batch.json")+",\"actual_parent_invocation\":"+quote(invocation_)+",\"manifest_sha256\":"+quote(manifestHash_)+",\"binary_sha256\":"+quote(binaryHash_)+",\"postgraph_flags\":"+flagJson(a.flags)+",\"effective_child_options\":"+config_+",\"report_total_seconds_scope\":\"postgraph_child_only; preprocessing reused from batch\",\"report_peak_memory_bytes_scope\":\"child_process_only; excludes parent preprocessing high-water\",\"postgraph_seconds\":"+number(seconds)+"}\n";
    return atomicWrite(root_+"/"+a.id+"/fork_child.json",text,error);
}
ForkBatch::Entry ForkBatch::enter(AssemblyOptions& o,const UnitigGraph& graph,const SequenceStore& reads,double prefixSeconds,std::string& error) {
    if(!enabled_)return Entry::Disabled;
    if(!reads.paired()){error="developer fork batch requires paired reads";return Entry::Failed;}
    prefixSeconds_=prefixSeconds;
    state_="{\"stage\":\"final_graph_pre_resolver\",\"k\":"+std::to_string(graph.k())+",\"graph_node_slots\":"+std::to_string(graph.nodes.size())+",\"live_nodes\":"+std::to_string(graph.liveCount())+",\"graph_bases\":"+std::to_string(graph.totalLength())+",\"reads\":"+std::to_string(reads.size())+",\"read_bases\":"+std::to_string(reads.totalBases())+",\"paired_reads\":"+std::to_string(reads.pairedReads())+",\"identity_basis\":\"same unchanged in-memory parent; no serialized checkpoint\"}";
    if(!writeBatch("running",error))return Entry::Failed;
    if (!safeSignalState(error)) { writeBatch("failed", error); return Entry::Failed; }
    SignalScope signals;interrupted=0;if(!signals.install()){error="cannot install batch interruption handling";return Entry::Failed;}
    bool anyFailure=false;
    for(size_t i=0;i<arms_.size();++i){
        auto& a=arms_[i];const std::string dir=root_+"/"+a.id;
        if(interrupted){error="developer batch interrupted";anyFailure=true;break;}
        if(taskCount()!=1){error="developer fork requires exactly one OS thread";anyFailure=true;break;}
        if(::mkdir(dir.c_str(),0755)!=0){error="cannot create isolated developer output: "+dir;anyFailure=true;break;}
        if(std::fflush(nullptr)!=0){error="cannot flush streams before developer fork";anyFailure=true;break;}
        auto start=std::chrono::steady_clock::now();pid_t pid=::fork();
        if(pid<0){error="developer fork failed: "+std::string(std::strerror(errno));anyFailure=true;break;}
        if(pid==0){
            signals.restore();
            // This standalone developer executable installs no child signal
            // handlers. Explicit defaults ensure parent cancellation/death can
            // terminate the child even if the launcher ignored these signals.
            struct sigaction childDefault{};
            childDefault.sa_handler = SIG_DFL;
            sigemptyset(&childDefault.sa_mask);
            if (sigaction(SIGINT, &childDefault, nullptr) != 0 ||
                sigaction(SIGTERM, &childDefault, nullptr) != 0) ::_exit(125);
            sigset_t childUnblock;
            ::sigemptyset(&childUnblock);
            ::sigaddset(&childUnblock, SIGINT);
            ::sigaddset(&childUnblock, SIGTERM);
            if (::sigprocmask(SIG_UNBLOCK, &childUnblock, nullptr) != 0) ::_exit(125);
            if(::prctl(PR_SET_PDEATHSIG,SIGTERM)!=0||::getppid()!=parentPid_)::_exit(125);
            int log=::open((dir+"/diagnostic.log").c_str(),O_CREAT|O_EXCL|O_WRONLY,0644);
            if(log<0||::dup2(log,STDOUT_FILENO)<0||::dup2(log,STDERR_FILENO)<0)::_exit(125);
            if(log>STDERR_FILENO)::close(log);
            for(const char* f:flags)if(::unsetenv(f)!=0)::_exit(125);
            for(const auto& f:a.flags)if(::setenv(f.first.c_str(),f.second.c_str(),1)!=0)::_exit(125);
            ::unsetenv("TESSERACT_DEV_FORK_BATCH");
            o.outDir=dir;child_=true;childIndex_=i;config_=optionJson(o);
            if(!writeChild("running",0,error))return Entry::Failed;
            std::fprintf(stderr,"[developer fork-child] arm=%s; shared preprocessing=%s/fork_batch.json; reported elapsed is postgraph only\n",a.id.c_str(),root_.c_str());
            return Entry::Child;
        }
        a.pid=pid;a.status="running";
        if(interrupted)::kill(pid,interrupted);
        int status=0;struct rusage usage{};bool waited=false;
        while(true){
            a.peakPss=std::max(a.peakPss,pss(parentPid_)+pss(pid));
            pid_t result=::wait4(pid,&status,WNOHANG,&usage);
            if(result==pid){waited=true;break;}
            if(result<0&&errno!=EINTR){
                const int waitError = errno;
                error="cannot reap developer child";
                // ECHILD means this PID is no longer known to be ours.
                if (waitError != ECHILD) {
                    ::kill(pid,SIGTERM);
                    while(::waitpid(pid,&status,0)<0&&errno==EINTR){}
                }
                break;
            }
            if(interrupted)::kill(pid,interrupted);
            struct timespec delay{0,20000000};::nanosleep(&delay,nullptr);
        }
        a.seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        a.userSeconds=usage.ru_utime.tv_sec+usage.ru_utime.tv_usec/1e6;a.systemSeconds=usage.ru_stime.tv_sec+usage.ru_stime.tv_usec/1e6;
        if(waited&&WIFEXITED(status))a.exitCode=WEXITSTATUS(status);
        if(waited&&WIFSIGNALED(status))a.signal=WTERMSIG(status);
        bool ok=waited&&a.exitCode==0&&exists(dir+"/contigs.fasta")&&exists(dir+"/report.json")&&exists(dir+"/fork_child.json");
        a.status=ok?"complete":"failed";anyFailure|=!ok;
        if(!writeBatch("running",error)){anyFailure=true;break;}
        if(interrupted){error="developer batch interrupted";anyFailure=true;break;}
    }
    if(!writeBatch(anyFailure?"failed":"complete",error))return Entry::Failed;
    if(anyFailure){if(error.empty())error="one or more developer fork arms failed; inspect fork_batch.json";return Entry::Failed;}
    std::fprintf(stderr,"[developer fork-batch complete] %zu diagnostic arms; no parent assembly; %s/fork_batch.json\n",arms_.size(),root_.c_str());
    return Entry::ParentComplete;
}
bool ForkBatch::finishChild(double seconds,std::string& error) {
    if(!child_)return true;
    if(!writeChild("complete",seconds,error))return false;
    childFinished_=true;return true;
}
ForkBatch::~ForkBatch(){if(child_&&!childFinished_){std::string error;if(!writeChild("failed",0,error))std::fprintf(stderr,"developer provenance failure: %s\n",error.c_str());}}
} } // namespace ts::dev
