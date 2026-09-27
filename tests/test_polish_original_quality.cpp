// Integration test of the default-off TESSERACT_POLISH_ORIGINAL_QUALITY polisher path.
//
// Ported from the untracked test_polish_original_quality.cpp.stale, which linked a
// symbol-renamed pre-change polish.cpp as ts::nativeBaselinePolish and so could not be
// built from the release tree. The native baseline is now polishContigs itself with the
// flag unset: polish.cpp reads the flag (quality_consensus::enabled()) at polish time, and
// the flag-unset path is exactly what the release ships. What this pins that nothing else
// does: every value other than a literal "1" is off and equals the native polisher, and a
// native correction keeps priority over a quality proposal that would oppose it.
#include "test_env.h"
#include "env_reject.h"
#include "polish.h"
#include "polish_quality.h"
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <random>
#include <string>
#include <unistd.h>
// Native baseline = polishContigs with TESSERACT_POLISH_ORIGINAL_QUALITY unset
// (polish.cpp gates the quality path on the flag at polish time), env restored after.
namespace { ts::PolishStats nativeBaselinePolish(std::vector<std::string>& c,const ts::SequenceStore& r,int t,int k,int d,double f){
    const char* v=std::getenv("TESSERACT_POLISH_ORIGINAL_QUALITY");const bool had=v!=nullptr;const std::string saved=had?v:"";
    unsetenv("TESSERACT_POLISH_ORIGINAL_QUALITY");auto s=ts::polishContigs(c,r,t,k,d,f);
    if(had)setenv("TESSERACT_POLISH_ORIGINAL_QUALITY",saved.c_str(),1);
    return s;} }
namespace {
size_t checks=0,serial=0;
std::filesystem::path dir;
void check(bool ok,const char* label){++checks;if(!ok){std::fprintf(stderr,"FAIL: %s\n",label);std::exit(1);}}
void flag(const char* value){if(value)setenv("TESSERACT_POLISH_ORIGINAL_QUALITY",value,1);else unsetenv("TESSERACT_POLISH_ORIGINAL_QUALITY");}
std::string rc(std::string s){std::reverse(s.begin(),s.end());for(char& c:s){int b=ts::baseCode(c);c=b<0?'N':ts::codeBase(3-b);}return s;}
std::string sequence(unsigned seed){std::mt19937 g(seed);std::string s(650,'A');for(char& c:s)c="ACGT"[g()%4];s[300]='C';return s;}
struct Read{std::string seq,qual;};
struct Pair{Read first,second;};
Read observation(const std::string& contig,char allele,int q,int strand,int start=225){
    const int center=300-start;
    Read r{contig.substr(start,151),std::string(151,'I')};r.seq[center]=allele;
    if(strand)r.seq=rc(r.seq);
    r.qual[strand?150-center:center]=char(q+33);return r;
}
void fastq(std::ostream& f,const Read& r,size_t id,int mate){f<<'@'<<id<<'/'<<mate<<'\n'<<r.seq<<"\n+\n"<<r.qual<<'\n';}
ts::SequenceStore load(const std::vector<Pair>& pairs,const char* value="1",int threads=2){
    flag(value);auto prefix=dir/std::to_string(serial++);std::string a=prefix.string()+".1.fq",b=prefix.string()+".2.fq";
    {std::ofstream x(a),y(b);for(size_t i=0;i<pairs.size();++i){fastq(x,pairs[i].first,i,1);fastq(y,pairs[i].second,i,2);}}
    ts::Library l;l.r1=a;l.r2=b;ts::SequenceStore s;std::string error;check(s.load({l},threads,error),"paired FASTQ load");return s;
}
std::vector<Pair> fixture(const std::string& s,int fragments=20,int highQ=30,int lowQ=2){
    std::vector<Pair> pairs;
    for(int i=0;i<fragments;++i){char allele=i<fragments/2?'A':'C';int q=allele=='A'?highQ:lowQ;pairs.push_back({observation(s,allele,q,0),observation(s,allele,q,1)});}return pairs;
}
void nativeStats(const ts::PolishStats& a,const ts::PolishStats& b){
    check(a.readsUsed==b.readsUsed && a.basesChanged==b.basesChanged && a.positionsCovered==b.positionsCovered &&
          a.lowCoveragePositions==b.lowCoveragePositions && a.meanDepth==b.meanDepth,"native counters identity");
}
}
int main(){
    testenv::clearTesseractEnv();  // first: flags are cached in statics on first use
    dir=std::filesystem::temp_directory_path()/("test-original-quality-"+std::to_string(getpid()));std::filesystem::create_directory(dir);
    const std::string s=sequence(91526);auto pairs=fixture(s);auto reads=load(pairs);
    check(reads.hasOriginalQualities() && reads.size()==40 && reads.pairedReads()==40,"original Q enabled with stable paired IDs");
    check(reads.originalQualityBytes()==reads.totalBases()+((reads.totalBases()+31)/32+(reads.totalBases()+63)/64)*8,"exact retained buffer bytes");
    check(reads.originalBaseAt(0,75)==0 && reads.originalBaseAt(1,75)==3 && reads.originalQualityAt(1,75)==30,"original alleles/Q remain in raw orientation");
    std::vector<std::string> before{s};auto ordinary=nativeBaselinePolish(before,reads,2,31,15,.9);
    check(before[0]==s && ordinary.basesChanged==0,"fixture has mixed unchanged native pile");
    for(const char* v:{static_cast<const char*>(nullptr),"0"}){
        auto off=load(pairs,v);check(!off.hasOriginalQualities() && off.originalQualityBytes()==0,"only literal1 retains input provenance");
        std::vector<std::string> actual{s};auto stats=ts::polishContigs(actual,off,3,31,15,.9);check(actual==before,"off output equals actual old native function");nativeStats(stats,ordinary);
    }
    // build_v3 (T16, envflags): a malformed value is no longer read as off; it is refused with
    // exit status 2 naming the flag (checked in a forked child), so "only literal 1 enables" holds.
    for(const char* v:{"true","01","1 ",""}){
        check(exitsWithEnvError("TESSERACT_POLISH_ORIGINAL_QUALITY",[&]{flag(v);(void)load(pairs,v);}),
              (std::string("malformed TESSERACT_POLISH_ORIGINAL_QUALITY='")+v+"' is refused (exit 2)").c_str());
        flag(nullptr);
    }
    flag("1");
    std::vector<std::string> enabled{s};auto on=ts::polishContigs(enabled,reads,3,31,15,.9);
    check(enabled[0][300]=='A' && on.basesChanged==1 && on.quality.extraChanges==1,"high-Q mixed evidence adds one proposal");
    check(on.quality.candidatePositions==1 && on.quality.maxFragmentDepth==20 && on.quality.eligibleObservations==40,"native site candidates and fragment depth correct");
    check(on.readsUsed==ordinary.readsUsed && on.positionsCovered==ordinary.positionsCovered && on.meanDepth==ordinary.meanDepth,"native mapping/coverage unchanged");
    for(int threads:{1,2,5}){std::vector<std::string> actual{s};auto st=ts::polishContigs(actual,reads,threads,31,15,.9);check(actual==enabled && st.quality.informativeFragmentSites==on.quality.informativeFragmentSites && st.quality.extraChanges==1,"thread-count deterministic consensus");}
    {auto few=load(fixture(s,10));std::vector<std::string> actual{s};auto st=ts::polishContigs(actual,few,2,31,15,.9);check(actual[0]==s && st.quality.maxFragmentDepth==10 && st.quality.belowDepth==1,"20 mate observations remain only10 informative fragments");}
    {std::vector<Pair> p;for(int i=0;i<20;++i)p.push_back({observation(s,'C',2,0),observation(s,'A',35,1)});auto x=load(p);std::vector<std::string> a{s};auto st=ts::polishContigs(a,x,2,31,15,.9);check(a[0]==s && st.quality.oneOrientation==1,"all alternate observations on one strand abstain");}
    {std::vector<Pair> p;for(int i=0;i<20;++i)p.push_back(i%2?Pair{observation(s,'A',35,0),observation(s,'C',2,1)}:Pair{observation(s,'C',2,0),observation(s,'A',35,1)});auto x=load(p);std::vector<std::string> a{s};auto st=ts::polishContigs(a,x,2,31,15,.9);check(a[0][300]=='A' && st.quality.maxFragmentDepth==20,"discordant opposite mates average once per fragment with both orientations");}
    {auto p=fixture(s,20,35,35);for(auto& pair:p){pair.first=observation(s,'A',35,0);pair.second=observation(s,'C',35,1);}auto x=load(p);std::vector<std::string> a{s};auto st=ts::polishContigs(a,x,2,31,15,.9);check(a[0]==s && st.quality.tied==1,"equally strong opposite alleles abstain");}
    {auto x=load(fixture(s,20,0,0));std::vector<std::string> a{s};auto st=ts::polishContigs(a,x,2,31,15,.9);check(a[0]==s && st.quality.uniformObservations==40 && st.quality.maxFragmentDepth==0,"Q0 mates do not create informative depth");}
    {auto x=load(fixture(s,20,93,2));std::vector<std::string> a{s};auto st=ts::polishContigs(a,x,2,31,15,.9);check(st.quality.cappedObservations==20 && a[0][300]=='A',"high rawQ retained then capped only for scoring");check(x.originalQualityAt(0,75)==93,"originalQ is not overwritten by cap");}
    {auto p=fixture(s);p.push_back({observation(s,'N',40,0),observation(s,'N',40,1)});auto x=load(p);x.maskRange(2,75,76);x.setBase(4,75,1);x.setBase(40,75,0);std::vector<std::string> a{s};auto st=ts::polishContigs(a,x,2,31,15,.9);check(st.quality.excludedOriginalAmbiguous==2 && st.quality.excludedMasked==1 && st.quality.excludedModified==1,"original N, current mask and net substitutions excluded separately");check(x.originalBaseAt(4,75)==0 && x.originalQualityAt(4,75)==30 && x.baseAt(4,75)==1,"EC writes never alter original provenance");check(x.originalBaseAt(40,75)==-1,"inputN never reinterpreted as packed A");}
    {std::vector<Pair> p;for(int i=0;i<20;++i)p.push_back({observation(s,i==0?'C':'A',i==0?40:1,0),observation(s,i==0?'C':'A',i==0?40:1,1)});auto x=load(p);std::vector<std::string> a{s},b{s};auto old=nativeBaselinePolish(a,x,2,31,15,.9);auto st=ts::polishContigs(b,x,3,31,15,.9);check(a==b && b[0][300]=='A' && old.basesChanged==1 && st.quality.extraChanges==0 && st.quality.candidatePositions==0,"native correction has priority even when quality score would oppose it");}
    {flag("1");auto file=dir/"noq.fa";std::ofstream f(file);for(size_t i=0;i<pairs.size();++i)f<<'>'<<i<<"a\n"<<pairs[i].first.seq<<'\n'<<'>'<<i<<"b\n"<<pairs[i].second.seq<<'\n';f.close();ts::Library l;l.r1=file.string();ts::SequenceStore x;std::string error;check(x.load({l},1,error),"FASTA provenance load");std::vector<std::string> a{s};auto st=ts::polishContigs(a,x,2,31,15,.9);check(a[0]==s && st.quality.excludedMissingQuality==40,"FASTA has no invented Phred confidence");}
    {auto p=fixture(s);auto x=load(p,nullptr);flag("1");std::vector<std::string> a{s};auto st=ts::polishContigs(a,x,2,31,15,.9);check(a[0]==s && !st.quality.provenanceAvailable && st.quality.extraChanges==0,"flag enabled after loader fails closed without original provenance");}
    {auto p=fixture(s);Read empty{std::string(100,'A'),std::string(100,'!')};p.insert(p.begin(),{empty,observation(s,'C',2,1)});auto x=load(p);check(x.length(0)==0 && x.length(1)==151 && x.originalBaseAt(2,75)==0,"fully trimmed first mate retains IDs");std::vector<std::string> a{s};auto st=ts::polishContigs(a,x,2,31,15,.9);check(a[0][300]=='A' && st.quality.maxFragmentDepth==21,"zero-length first mate preserves fragment grouping");}
    {auto p=fixture(s);for(auto& pair:p){pair.second.qual.replace(147,4,std::string(4,char(20+33)));pair.second.seq+=std::string(20,'G');pair.second.qual+=std::string(20,'!');}auto x=load(p);check(x.length(1)==151 && x.trimmedBases()==400,"reverse reads have asymmetric raw and retained lengths");std::vector<std::string> a{s};auto st=ts::polishContigs(a,x,2,31,15,.9);check(a[0][300]=='A' && st.quality.extraChanges==1 && st.quality.excludedModified==0,"reverse coordinate uses retained length");}
    {std::vector<std::string> a{s,s};auto st=ts::polishContigs(a,reads,2,31,15,.9);check(a[0]==s && a[1]==s && st.readsUsed==0 && st.quality.candidatePositions==0,"duplicate native copies do not gain pooled evidence");}
    {const auto other=sequence(74151);std::vector<Pair> p;for(int i=0;i<20;++i){char allele=i<10?'A':'C';int q=allele=='A'?30:2;p.push_back({observation(s,allele,q,i%2),observation(other,allele,q,i%2)});}auto x=load(p);std::vector<std::string> a{s,other};auto st=ts::polishContigs(a,x,3,31,15,.9);check(a[0][300]=='A' && a[1][300]=='A' && st.quality.extraChanges==2 && st.quality.maxFragmentDepth==20 && st.quality.informativeFragmentSites==40,"mates on distinct contig copies scored separately without pooled depth");}
    for (const auto& crop : {std::pair<int,int>{280,140},std::pair<int,int>{180,141}}) {
        std::vector<Pair> p;for(int i=0;i<20;++i){char allele=i<10?'A':'C';int q=allele=='A'?30:2;p.push_back({observation(s,allele,q,0,230),observation(s,allele,q,1,230)});}
        auto x=load(p);std::vector<std::string> a{s.substr(crop.first,crop.second)};
        auto st=ts::polishContigs(a,x,3,31,15,.9);
        check(a[0][300-crop.first]=='A' && st.quality.eligibleObservations==40 && st.quality.excludedModified==0,"left/right overhang clipping preserves asymmetric forward/reverse original coordinates");
    }
    // Paired-first mixed-library layout keeps quality offsets attached to input records.
    {auto paired=load(pairs);auto single=dir/"single.fq";{std::ofstream f(single);fastq(f,observation(s,'A',27,0),1,1);}auto a=dir/(std::to_string(serial-1)+".1.fq"),b=dir/(std::to_string(serial-1)+".2.fq");ts::Library se,pe;se.r1=single.string();pe.r1=a.string();pe.r2=b.string();ts::SequenceStore x;std::string error;check(x.load({se,pe},3,error),"mixed single/paired load");check(x.pairedReads()==40 && !x.hasMate(40) && x.originalQualityAt(40,75)==27 && x.originalQualityAt(0,75)==30,"paired-first ordering preserves original qualities");flag("0");check(x.load({se},1,error) && !x.hasOriginalQualities() && x.originalQualityBytes()==0,"reloading disabled releases provenance buffers");flag("1");se.r1=(dir/"missing.fastq").string();check(!x.load({se},1,error) && !x.hasOriginalQualities(),"load failure cannot retain stale provenance");}
    std::filesystem::remove_all(dir);
    std::printf("%zu original-quality integration checks passed\n",checks);
}
