// Actual graph anchoring, internal insert fitting and reciprocal resolution.
#include "resolve.h"
#include <algorithm>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <unistd.h>

namespace {
std::string dna(size_t n, unsigned seed) {
    std::mt19937 rng(seed); std::string s(n,'A');
    for (char& c : s) c = "ACGT"[rng()%4];
    return s;
}
void edge(ts::UnitigGraph& g, unsigned a, unsigned b) {
    g.nodes[a].ends[1].push_back({b,0}); g.nodes[b].ends[0].push_back({a,1});
}
struct Temp {
    std::filesystem::path path;
    explicit Temp(int id) : path(std::filesystem::temp_directory_path()/
        ("tesseract_route_test_"+std::to_string(getpid())+"_"+std::to_string(id))) {
        std::filesystem::create_directory(path);
    }
    ~Temp() { std::filesystem::remove_all(path); }
};
struct Result { int arm = -1; size_t length = 0; std::vector<std::string> contigs; };
}
int main() {
    setenv("TESSERACT_COMMON_PREFIX","0",1);
    setenv("TESSERACT_JOIN_TRACE","1",1);
    for (auto flag : {"TESSERACT_WEIGHTED_RESOLVER_COVERAGE", "TESSERACT_WEIGHTED_ELIGIBLE_COVERAGE",
                      "TESSERACT_PREFIX_SNP_BUBBLES", "TESSERACT_EXCLUDE_SHARED_REPEAT_SUPPORT",
                      "TESSERACT_SHARED_SUPPORT_AUDIT"}) unsetenv(flag);
    int checks = 0, cases = 0;
    auto check = [&](bool good, const char* why) { ++checks; if (!good) throw std::runtime_error(why); };
    auto scenario = [&](bool correctLong, bool sameLength, bool reorder, bool reverse, int training, bool forcedWindow = true, int comparablePairs = 41) {
        const int id = ++cases;
        ts::UnitigGraph original; original.setK(31); original.nodes.resize(9);
        auto& n = original.nodes;
        const int endpointLength = comparablePairs < 41 ? 5000 : 1000;
        n[0].seq=dna(endpointLength,100);
        const std::string suffix=dna(30,101);
        n[1].seq=n[0].seq.substr(endpointLength-30)+dna(70,102)+suffix;
        n[2].seq=n[0].seq.substr(endpointLength-30)+dna(sameLength?70:170,103)+suffix;
        n[3].seq=suffix+dna(endpointLength-30,104); n[4].seq=dna(4000,105);
        for (size_t i=5;i<n.size();++i) n[i].seq=dna(200,200+i);
        for (auto& node:n) node.coverage=10;
        n[1].coverage=correctLong?60:40;
        n[2].coverage=correctLong?40:60;
        edge(original,0,1); edge(original,1,3); edge(original,0,2); edge(original,2,3);
        const std::string shortPath=n[0].seq+n[1].seq.substr(30)+n[3].seq.substr(30);
        const std::string longPath=n[0].seq+n[2].seq.substr(30)+n[3].seq.substr(30);
        Temp temp(id);
        const auto r1=temp.path/"r1.fastq",r2=temp.path/"r2.fastq";
        std::ofstream f(r1),r(r2); size_t count=0, trained=0;
        auto pair=[&](unsigned a,int pos,unsigned b,int pos2) {
            const auto s1=n[a].seq.substr(pos,70),s2=ts::reverseComplement(n[b].seq.substr(pos2,70));
            const auto name="@p"+std::to_string(count++);
            f<<name<<"/1\n"<<s1<<"\n+\n"<<std::string(70,'I')<<'\n';
            r<<name<<"/2\n"<<s2<<"\n+\n"<<std::string(70,'I')<<'\n';
        };
        if (training) for (int z=-10;z<=10;++z) {
            const int fragment=(correctLong?600:500)+z*10;
            for (int rep=0;rep<20*(11-std::abs(z));++rep) {
                if (training>0 && int(trained)>=training) break;
                const int start=100+((trained*37)%2000);
                pair(4,start,4,start+fragment-70); ++trained;
            }
        }
        for (int i=0;i<comparablePairs;++i) pair(0,endpointLength-200+i-comparablePairs/2,3,160);
        // Valid long-fragment endpoint links inside the forced legacy window
        // but longer than every calibration interval; they cannot certify a
        // route density. Keep total endpoint support at 41 in sparse cases.
        for (int i=comparablePairs;i<41;++i) pair(0,500+i,3,3500);
        f.close();r.close();
        ts::Library lib;lib.r1=r1.string();lib.r2=r2.string();
        ts::SequenceStore reads; std::string error;
        if (!reads.load({lib},1,error)) throw std::runtime_error(error);
        ts::UnitigGraph g=original;
        if (reorder) {
            const std::vector<unsigned> perm{8,2,6,1,3,7,5,4,0};
            for (size_t i=0;i<n.size();++i) {
                g.nodes[perm[i]]=n[i];
                for (auto& links:g.nodes[perm[i]].ends) for (auto& link:links) link.to=perm[link.to];
            }
        }
        if (reverse) for (auto& node:g.nodes) {
            node.seq=ts::reverseComplement(node.seq); std::swap(node.ends[0],node.ends[1]);
            for (auto& links:node.ends) for (auto& link:links) link.toEnd=1-link.toEnd;
        }
        check(g.validate().empty(),"invalid transformed graph");
        auto run=[&](int flag) {
            if (flag<0) unsetenv("TESSERACT_ROUTE_DISTANCE");
            else setenv("TESSERACT_ROUTE_DISTANCE",flag?"1":"0",1);
            ts::PairedResolver resolver(g,reads,1,2,1.02,.02);
            if (forcedWindow && training>=0 && training<1000) resolver.setInsertBounds(0,1000);
            if (comparablePairs < 41) resolver.setInsertBounds(0,10000);
            resolver.buildSupport();
            check(resolver.stats().pairsLinking==41,"spanning pairs lost in anchoring");
            std::vector<double> cov;Result result;resolver.resolve(result.contigs,cov);
            for (const auto& s:result.contigs) {
                if (s==shortPath || ts::reverseComplement(s)==shortPath) { result.arm=0;result.length=s.size(); }
                if (s==longPath || ts::reverseComplement(s)==longPath) { result.arm=1;result.length=s.size(); }
            }
            check(result.arm>=0,"unique flanks no longer joined through one arm");
            return result;
        };
        const auto baseline=run(-1),explicitOff=run(0),candidate=run(1);
        check(baseline.contigs==explicitOff.contigs,"default-off identity failed");
        check(baseline.arm==(correctLong?0:1),"fixture failed to pick high-coverage baseline arm");
        check(baseline.contigs.size()==candidate.contigs.size(),"route model changes contig count");
        if (sameLength || (training>=0 && training<1000) || comparablePairs<=2)
            check(baseline.contigs==candidate.contigs,"unidentifiable route did not retain exact legacy output");
        else {
            check(candidate.arm==(correctLong?1:0),"distance evidence chose wrong route");
            check(candidate.length==size_t(2*endpointLength+(correctLong?200:100)-30),"wrong assembled path length");
        }
    };
    scenario(true,false,false,false,-1); // true long connector, low repeat coverage
    scenario(false,false,false,false,-1); // true short connector, low repeat coverage
    scenario(true,false,true,false,-1); // node-ID permutation
    scenario(true,false,true,true,-1); // every graph node reverse complemented
    scenario(false,false,true,true,-1); // reciprocal short-correct evidence
    scenario(false,true,false,false,-1); // sequence alternatives have identical length
    scenario(true,false,false,false,0); // forced broad window but no internal training
    scenario(true,false,false,false,100); // insufficient internal training
    scenario(true,false,false,false,0,false); // no usable insert model at all
    scenario(true,false,false,false,-1,true,1); // one informative physical pair cannot change route
    scenario(true,false,false,false,-1,true,2); // two fractional contrasts remain below bar2
    scenario(true,false,false,false,-1,true,3); // three strong contrasts clear the existing bar
    // After A-X and Y-B merge, the forward history A->Y favors the long
    // interior while reciprocal B->X favors the short one. Endpoint support
    // remains mutual; the model must restore legacy paths without refusing it.
    auto mirrorScenario = [&](bool reorder, bool reverse) {
        const int id=++cases;
        ts::UnitigGraph original;original.setK(31);original.nodes.resize(12);
        auto& n=original.nodes;
        n[0].seq=dna(1000,300); // A
        n[1].seq=n[0].seq.substr(970)+dna(40,301); // X, 70 bases
        const std::string suffix=dna(30,302);
        n[2].seq=n[1].seq.substr(40)+dna(70,303)+suffix;
        n[3].seq=n[1].seq.substr(40)+dna(170,304)+suffix;
        n[4].seq=suffix+dna(40,305); // Y
        n[5].seq=n[4].seq.substr(40)+dna(970,306); // B
        n[6].seq=dna(4000,307);
        for(size_t i=7;i<n.size();++i)n[i].seq=dna(200,400+i);
        for(auto& node:n)node.coverage=10;
        n[2].coverage=60;n[3].coverage=40;
        edge(original,0,1);edge(original,1,2);edge(original,1,3);
        edge(original,2,4);edge(original,3,4);edge(original,4,5);
        const std::string expected=n[0].seq+n[1].seq.substr(30)+n[2].seq.substr(30)+
                                   n[4].seq.substr(30)+n[5].seq.substr(30);
        Temp temp(id);const auto r1=temp.path/"r1.fastq",r2=temp.path/"r2.fastq";
        std::ofstream f(r1),r(r2);size_t count=0,trained=0;
        auto pair=[&](unsigned a,int p1,unsigned b,int p2) {
            const auto s1=n[a].seq.substr(p1,70),s2=ts::reverseComplement(n[b].seq.substr(p2,70));
            const auto name="@m"+std::to_string(count++);
            f<<name<<"/1\n"<<s1<<"\n+\n"<<std::string(70,'I')<<'\n';
            r<<name<<"/2\n"<<s2<<"\n+\n"<<std::string(70,'I')<<'\n';
        };
        for(int z=-10;z<=10;++z)for(int rep=0;rep<20*(11-std::abs(z));++rep) {
            const int start=100+((trained*37)%2000);
            pair(6,start,6,start+600+10*z-70);++trained;
        }
        for(int i=0;i<100;++i){pair(0,500+i%5,1,0);pair(4,0,5,400+i%5);}
        for(int p1=660;p1<=700;++p1)pair(0,p1,4,0); // long favored after X history shift
        for(int p2=330;p2<=370;++p2)pair(1,0,5,p2); // short favored after Y history shift
        f.close();r.close();
        ts::Library lib;lib.r1=r1.string();lib.r2=r2.string();ts::SequenceStore reads;std::string error;
        if(!reads.load({lib},1,error))throw std::runtime_error(error);
        ts::UnitigGraph g=original;
        if(reorder) {
            const std::vector<unsigned> perm{11,7,4,2,8,0,3,1,5,10,6,9};
            for(size_t i=0;i<n.size();++i){g.nodes[perm[i]]=n[i];
                for(auto& links:g.nodes[perm[i]].ends)for(auto& link:links)link.to=perm[link.to];}
        }
        if(reverse)for(auto& node:g.nodes){node.seq=ts::reverseComplement(node.seq);
            std::swap(node.ends[0],node.ends[1]);
            for(auto& links:node.ends)for(auto& link:links)link.toEnd=1-link.toEnd;}
        check(g.validate().empty(),"invalid mirror-history graph");
        auto run=[&](bool enabled){setenv("TESSERACT_ROUTE_DISTANCE",enabled?"1":"0",1);
            ts::PairedResolver resolver(g,reads,1,2,1.02,.02);resolver.buildSupport();
            check(resolver.stats().pairsLinking==282,"mirror fixture loses linked pairs");
            std::vector<std::string> result;std::vector<double> cov;resolver.resolve(result,cov);return result;};
        const auto baseline=run(false),candidate=run(true);
        check(baseline==candidate,"contradictory mirror histories do not restore exact legacy output");
        bool joined=false;
        for(const auto& s:candidate)if(s==expected||ts::reverseComplement(s)==expected)joined=true;
        check(joined,"mirror disagreement rejected or shortened the legacy join");
    };
    mirrorScenario(false,false);
    mirrorScenario(true,true);
    std::cout<<checks<<" full resolver route checks passed across "<<cases<<" fixtures\n";
}
