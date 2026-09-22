// Standalone focused test. Compile with correct.cpp, seqio.cpp, counter.cpp, util.cpp.
#include "correct.h"
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <random>
#include <string>
#include <unistd.h>

namespace {
int checks=0;
void check(bool ok,const char* label) { ++checks; if(!ok){std::fprintf(stderr,"FAIL: %s\n",label);std::exit(1);} }
std::string rc(std::string s) {
    std::reverse(s.begin(),s.end());
    for(char& c:s) c=c=='N'?'N':ts::codeBase(3-ts::baseCode(c));
    return s;
}
std::filesystem::path directory;
std::mt19937 generator(20260919);
std::string randomSeq(int n) {std::string s(n,'A');for(char& c:s)c="ACGT"[generator()%4];return s;}
void add(ts::KmerTable& solid,const std::string& s,int k,unsigned depth,int first=0,int last=100000) {
    ts::Kmer f;
    for(size_t i=0;i<s.size();++i) {
        f=ts::pushBack(f,ts::baseCode(s[i]),k);
        const int start=int(i)-k+1;
        if(start>=first && start<=last) solid.put(ts::canonical(f,k),depth);
    }
}
ts::SequenceStore load(const std::string& s) {
    std::ofstream(directory/"reads.fa")<<">forward\n"<<s<<"\n>reverse\n"<<rc(s)<<'\n';
    ts::Library library;library.r1=(directory/"reads.fa").string();
    ts::SequenceStore reads;std::string error;
    check(reads.load({library},1,error),"load fixture");return reads;
}
void flag(const char* v) {if(v)setenv("TESSERACT_EC_REQUIRE_UNIQUE_BEST",v,1);else unsetenv("TESSERACT_EC_REQUIRE_UNIQUE_BEST");}
void same(const ts::SequenceStore& a,const ts::SequenceStore& b,const char* label) {
    check(a.size()==b.size(),label);
    for(size_t i=0;i<a.size();++i)check(a.decode(i)==b.decode(i),label);
}
void sameStats(const ts::CorrectionStats& a,const ts::CorrectionStats& b) {
    check(a.readsExamined==b.readsExamined && a.readsCorrected==b.readsCorrected &&
          a.basesCorrected==b.basesCorrected && a.readsUncorrectable==b.readsUncorrectable &&
          a.basesMasked==b.basesMasked,"legacy counters unchanged");
}
}
int main() {
    directory=std::filesystem::path("/tmp")/("test-ec-unique-"+std::to_string(getpid()));
    std::filesystem::create_directory(directory);
    for(int k:{21,33})for(int position:{25,90})for(bool unequal:{false,true}) {
        std::string raw=randomSeq(150);raw[position]='C';
        std::string g=raw,t=raw;g[position]='G';t[position]='T';
        ts::KmerTable solid;add(solid,g,k,2);add(solid,t,k,unequal?50:2);
        auto legacy=load(raw);flag(nullptr);const auto old=ts::correctReads(legacy,solid,k,1,1000);
        check(legacy.decode(0)==g && rc(legacy.decode(1))==t,"off retains exact alphabet-order behavior");
        check(old.basesCorrected==2 && old.ambiguousExtensions==0 && old.ambiguityMaskedBases==0,"off counters");
        for(const char* v:{"","0","true","01","1x"}) {
            auto disabled=load(raw);flag(v);const auto ds=ts::correctReads(disabled,solid,k,2,1000);
            same(disabled,legacy,"only exact literal 1 enables");sameStats(ds,old);
            check(ds.ambiguousExtensions==0 && ds.ambiguityMaskedBases==0,"disabled diagnostic counters zero");
        }
        auto retained=load(raw);flag("1");const auto rs=ts::correctReads(retained,solid,k,2,1000);
        check(retained.decode(0)==raw && rc(retained.decode(1))==raw,"equal best abstains in both orientations");
        check(rs.basesCorrected==0 && rs.ambiguousExtensions==2 && rs.basesMasked==0 && rs.ambiguityMaskedBases==0,
              "above-length mask threshold retains original bases");
        const int maskLength=position==25?position+1:150-position;
        auto masked=load(raw);std::vector<uint32_t> maskIds;
        const auto ms=ts::correctReads(masked,solid,k,2,maskLength,&maskIds);
        std::string expected=raw;
        if(position==25)expected.replace(0,maskLength,maskLength,'N');else expected.replace(position,maskLength,maskLength,'N');
        check(masked.decode(0)==expected && rc(masked.decode(1))==expected,"existing inclusive minMaskRun boundary and orientation");
        check(ms.basesCorrected==0 && ms.ambiguousExtensions==2 && ms.basesMasked==size_t(2*maskLength) &&
              ms.ambiguityMaskedBases==ms.basesMasked,"mask diagnostics count only ambiguous stops");
        check(maskIds.size()==2,"mask provenance retained");
        for(size_t i=0;i<raw.size();++i)check(masked.rawBaseAt(0,i)==ts::baseCode(raw[i]),"mask preserves raw bases");
        auto boundary=load(raw);const auto bs=ts::correctReads(boundary,solid,k,1,maskLength+1);
        check(boundary.decode(0)==raw && rc(boundary.decode(1))==raw && bs.basesMasked==0,"below minimum tail unchanged");
    }
    // A weaker alternate is considered first alphabetically. The unique longer
    // corroborating run must still win, without frequency weighting.
    for(int position:{25,90}) {
        const int k=21;std::string raw=randomSeq(150);raw[position]='C';
        std::string g=raw,a=raw;g[position]='G';a[position]='A';
        ts::KmerTable solid;add(solid,g,k,2);
        if(position==25)add(solid,a,k,50,position-2,position);
        else add(solid,a,k,50,position-k+1,position-k+3);
        auto old=load(raw),now=load(raw);flag(nullptr);const auto os=ts::correctReads(old,solid,k,1);
        flag("1");const auto ns=ts::correctReads(now,solid,k,2);
        check(now.decode(0)==g && rc(now.decode(1))==g,"unique best unchanged even against deeper weaker candidate");
        same(old,now,"unique best off/on identity");sameStats(os,ns);
        check(ns.ambiguousExtensions==0 && ns.ambiguityMaskedBases==0,"unique best no ambiguity diagnostics");
    }
    for(int position:{1,148}) {
        const int k=21;std::string raw=randomSeq(150);raw[position]='C';
        std::string g=raw,t=raw;g[position]='G';t[position]='T';
        ts::KmerTable solid;add(solid,g,k,2);add(solid,t,k,2);
        auto reads=load(raw);flag("1");const auto st=ts::correctReads(reads,solid,k,2,8);
        check(reads.decode(0)==raw && rc(reads.decode(1))==raw,"short terminal ambiguity retains original below mask minimum");
        check(st.ambiguousExtensions==2 && st.basesCorrected==0 && st.basesMasked==0,"terminal corroboration tie is rejected");
    }
    {
        const int k=21,position=90;std::string raw=randomSeq(150);raw[position]='C';
        ts::KmerTable solid;
        for(char base:std::string("AGT")){std::string alt=raw;alt[position]=base;add(solid,alt,k,2);}
        auto reads=load(raw);flag("1");const auto st=ts::correctReads(reads,solid,k,2,1000);
        check(reads.decode(0)==raw && rc(reads.decode(1))==raw,"all three alternatives tied abstain");
        check(st.ambiguousExtensions==2 && st.basesCorrected==0,"three tied candidates count one event per extension");
    }
    for(int position:{25,90}) {
        const int k=21;std::string raw=randomSeq(150);raw[position]='T';
        std::string g=raw;g[position]='G';ts::KmerTable solid;add(solid,g,k,2);
        for(char base:std::string("AC")) {
            std::string alt=raw;alt[position]=base;
            if(position==25)add(solid,alt,k,50,position-2,position);
            else add(solid,alt,k,50,position-k+1,position-k+3);
        }
        auto reads=load(raw);flag("1");const auto st=ts::correctReads(reads,solid,k,2);
        check(reads.decode(0)==g && rc(reads.decode(1))==g,"strictly better candidate resets an earlier lower-score tie");
        check(st.basesCorrected==2 && st.ambiguousExtensions==0,"lower-score tie does not reject unique maximum");
    }
    {
        const int k=21,position=90;std::string raw=randomSeq(150);raw[position]='T';
        ts::KmerTable solid;add(solid,raw,k,2,0,position-k);add(solid,raw,k,2,position+1);
        for(char base:std::string("AC")){std::string alt=raw;alt[position]=base;add(solid,alt,k,2,position-k+1,position-k+3);}
        auto old=load(raw),now=load(raw);flag(nullptr);const auto os=ts::correctReads(old,solid,k,1);
        flag("1");const auto ns=ts::correctReads(now,solid,k,2);
        same(old,now,"ties below corroboration requirement retain ordinary failure");sameStats(os,ns);
        check(ns.ambiguousExtensions==0 && ns.ambiguityMaskedBases==0 && ns.basesMasked>0,"inadequate tie is not an ambiguity event");
    }
    {
        const int k=21;std::string raw=randomSeq(150);raw[25]=raw[120]='C';
        ts::KmerTable solid;
        for(char left:std::string("GT"))for(char right:std::string("GT")) {
            std::string alt=raw;alt[25]=left;alt[120]=right;add(solid,alt,k,2);
        }
        auto reads=load(raw);flag("1");std::vector<uint32_t> ids;
        const auto st=ts::correctReads(reads,solid,k,2,8,&ids);
        std::string expected=raw;expected.replace(0,26,26,'N');expected.replace(120,30,30,'N');
        check(reads.decode(0)==expected && rc(reads.decode(1))==expected,"two sides of one trusted anchor have disjoint masks");
        check(st.ambiguousExtensions==4 && st.basesMasked==112 && st.ambiguityMaskedBases==112,"two events per physical read count exact masked union");
        check(ids.size()==4,"two mask provenance records per physical read retain existing convention");
        raw[50]='N';auto ambiguous=load(raw);const auto ns=ts::correctReads(ambiguous,solid,k,2);
        check(ambiguous.decode(0)==raw && rc(ambiguous.decode(1))==raw,"input N remains uncorrected under existing skip policy");
        check(ns.ambiguousExtensions==0 && ns.basesCorrected==0 && ns.basesMasked==0,"input N does not create tie diagnostics");
        raw[50]='A';auto masked=load(raw);masked.maskRange(0,50,51);masked.maskRange(1,99,100);
        const auto before0=masked.decode(0),before1=masked.decode(1);const auto ms=ts::correctReads(masked,solid,k,2);
        check(masked.decode(0)==before0 && masked.decode(1)==before1,"prior mask preserves existing skip policy");
        check(ms.ambiguousExtensions==0 && ms.basesCorrected==0 && ms.basesMasked==0,"prior mask does not create tie diagnostics");
    }
    flag(nullptr);std::filesystem::remove_all(directory);
    std::printf("%d correction unique-best checks passed\n",checks);
}
