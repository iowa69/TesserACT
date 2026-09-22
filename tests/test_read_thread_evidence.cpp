// Reference-free tests of exact routes and physical-molecule accounting.
#include <algorithm>
#include <cstdio>
#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <random>
#include <stdexcept>
#include <string>
#include <unistd.h>
#include <vector>
#include "read_thread_evidence.h"
#include "graph.h"
#include "seqio.h"
namespace {
int checks=0;
void check(bool ok,const char* why){++checks;if(!ok)throw std::runtime_error(why);}
std::string dna(size_t n,unsigned seed){std::mt19937 r(seed);std::string s(n,'A');for(char& c:s)c="ACGT"[r()%4];return s;}
void edge(ts::UnitigGraph& g,uint32_t a,uint32_t b){g.nodes[a].ends[1].push_back({b,0});g.nodes[b].ends[0].push_back({a,1});}
struct Fixture {
 ts::UnitigGraph graph;std::vector<uint8_t> eligible{1,1,0,1,1};
 std::string left=dna(500,11),otherLeft=dna(500,12),repeat=dna(140,13),right=dna(500,14),otherRight=dna(500,15);
 explicit Fixture(bool identical=false,int k=127){
  graph.setK(k);if(identical)otherRight=right;size_t ov=size_t(k-1);
  for(const auto& s:{left+repeat.substr(0,ov),otherLeft+repeat.substr(0,ov),repeat,repeat.substr(repeat.size()-ov)+right,repeat.substr(repeat.size()-ov)+otherRight}){
   ts::Unitig n;n.seq=s;n.coverage=graph.nodes.size()==2?40:20;graph.nodes.push_back(n);
  }
  edge(graph,0,2);edge(graph,1,2);edge(graph,2,3);edge(graph,2,4);check(graph.validate().empty(),"bidirected fixture");
 }
 std::pair<std::string,std::string> pair(int n,bool alt=false)const{
  const std::string f=left.substr(left.size()-size_t(n-156))+repeat+(alt?otherRight:right).substr(0,36);
  return {f.substr(0,n),ts::reverseComplement(f.substr(20,n))};
 }
};
struct Temp{std::filesystem::path dir=std::filesystem::temp_directory_path()/("tesseract-thread-"+std::to_string(getpid()));Temp(){std::filesystem::create_directories(dir);}~Temp(){std::filesystem::remove_all(dir);}};
ts::SequenceStore reads(const std::vector<std::pair<std::string,std::string>>& paired,const std::vector<std::string>& single={}){
 Temp temp;std::vector<ts::Library> libs;
 if(!paired.empty()){
  auto p1=temp.dir/"r1.fa",p2=temp.dir/"r2.fa";std::ofstream a(p1),b(p2);
  for(size_t i=0;i<paired.size();++i){a<<'>'<<i<<"/1\n"<<paired[i].first<<'\n';b<<'>'<<i<<"/2\n"<<paired[i].second<<'\n';}
  ts::Library l;l.r1=p1.string();l.r2=p2.string();libs.push_back(l);
 }
 if(!single.empty()){
  auto p=temp.dir/"single.fa";std::ofstream a(p);for(size_t i=0;i<single.size();++i)a<<'>'<<i<<'\n'<<single[i]<<'\n';
  ts::Library l;l.r1=p.string();libs.push_back(l);
 }
 ts::SequenceStore out;std::string error;bool loaded=out.load(libs,1,error);check(loaded,error.c_str());return out;
}
ts::ReadThreadPath canonical(ts::ReadThreadPath p){return std::min(p,ts::reverseReadThreadPath(p));}
auto collect(const Fixture& f,const ts::SequenceStore& r,ts::ReadThreadLimits lim={}){return ts::collectReadThreadEvidence(f.graph,r,f.eligible,lim);}
bool one(const ts::ReadThreadEvidence& e,ts::ReadThreadPath p,size_t count=1){return e.routes.size()==1&&e.routes[0].oriented==canonical(p)&&e.routes[0].fragments.size()==count;}
}
int main(int argc,char** argv){
 if(argc==3&&std::string(argv[1])=="--benchmark-pairs"){
  const size_t n=size_t(std::strtoull(argv[2],nullptr,10));check(n>0&&n<=100000,"bounded benchmark size");
  Fixture fixture;auto input=reads(std::vector<std::pair<std::string,std::string>>(n,fixture.pair(301)));
  const auto start=std::chrono::steady_clock::now();auto evidence=collect(fixture,input);
  const double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
  check(one(evidence,{0,4,6},n),"bounded throughput count");
  std::printf("{\"pairs\":%zu,\"reads\":%zu,\"bases\":%zu,\"seconds\":%.6f,\"fragments\":%zu}\n",n,2*n,602*n,seconds,evidence.stats.fragmentsAccepted);
  return 0;
 }
 for(int len:{251,301}){
  Fixture f;auto p=f.pair(len);auto e=collect(f,reads({p}));
  check(one(e,{0,4,6}),"exact source-interior-target route");
  check(e.stats.readsAccepted==2&&e.stats.fragmentsAccepted==1&&e.stats.matesDeduplicated==1,"overlapping mates count once");
  check(e.routes[0].fragments==std::vector<size_t>{0},"stable physical ID");
  auto repeat=collect(f,reads({p,p,p}));check(one(repeat,{0,4,6},3)&&repeat.routes[0].fragments==std::vector<size_t>({0,1,2}),"three molecules count three, not six or one");
  check(one(collect(f,reads({}, {p.first})),{0,4,6}),"single-end phase without insert model");
  auto unknown=collect(f,reads({{p.first,std::string(size_t(len),'N')}}));check(one(unknown,{0,4,6})&&unknown.stats.readsUnknown==1,"valid mate survives unavailable mate");
  check(one(collect(f,reads({}, {ts::reverseComplement(p.first)})),{0,4,6}),"reverse-complement read invariant");
  auto masked=reads({}, {p.first});masked.maskRange(0,125,126);auto me=collect(f,masked);check(me.routes.empty()&&me.stats.readsUnknown==1,"masked bases never converted to raw A");
  std::string altered=p.first;size_t m=size_t(len-156)+70;altered[m]=altered[m]=='A'?'C':'A';auto wrong=collect(f,reads({}, {altered}));check(wrong.routes.empty()&&wrong.stats.readsNoExactPath==1,"exact connector mismatch refused");
  Fixture amb(true);auto ae=collect(amb,reads({}, {p.first}));check(ae.routes.empty()&&ae.stats.readsAmbiguous==1,"two exact full routes refused");
  check(collect(f,reads({}, {p.first.substr(0,size_t(len-36))})).routes.empty(),"unspanned repeat is not phasing");
  auto alternative=f.pair(len,true);auto conflict=collect(f,reads({{p.first,alternative.second}}));check(conflict.routes.empty()&&conflict.stats.fragmentsConflicting==1,"conflicting mate routes refused");
  auto competing=collect(f,reads({p,alternative}));check(competing.routes.size()==2&&competing.stats.fragmentsAccepted==2,"alternative molecules remain visible");
  auto mix=collect(f,reads({p},{p.first}));check(one(mix,{0,4,6},2)&&mix.routes[0].fragments==std::vector<size_t>({0,1}),"single and pair IDs cannot collide");
  auto limited=collect(f,reads({}, {p.first}),{3,64});check(limited.routes.empty()&&limited.stats.readsSearchLimited==1,"exhaustion after first successful route fails closed");
  auto shallow=collect(f,reads({}, {p.first}),{4096,2});check(shallow.routes.empty()&&shallow.stats.readsSearchLimited==1,"node cap fails closed");
  Fixture overlap=f;overlap.graph.nodes[3].seq[0]=overlap.graph.nodes[3].seq[0]=='A'?'C':'A';check(collect(overlap,reads({}, {p.first})).routes.empty(),"invalid overlap refused");
  Fixture reciprocal=f;reciprocal.graph.nodes[3].ends[0].clear();check(collect(reciprocal,reads({}, {p.first})).routes.empty(),"nonreciprocal link refused");
  Fixture deleted=f;deleted.graph.nodes[3].deleted=true;check(collect(deleted,reads({}, {p.first})).routes.empty(),"deleted endpoint refused");
  Fixture disallowed=f;disallowed.eligible[3]=0;check(collect(disallowed,reads({}, {p.first})).routes.empty(),"caller endpoint eligibility respected");
  // Relabel every node and reverse selected stored unitig sequences.
  std::vector<uint32_t> ids{4,3,2,1,0};std::vector<uint8_t> flip{1,0,1,0,1},eligible(5);
  ts::UnitigGraph tg;tg.setK(127);tg.nodes.resize(5);
  for(size_t i=0;i<5;++i){
   auto& n=tg.nodes[ids[i]];n.seq=flip[i]?ts::reverseComplement(f.graph.nodes[i].seq):f.graph.nodes[i].seq;n.coverage=f.graph.nodes[i].coverage;eligible[ids[i]]=f.eligible[i];
   for(int end=0;end<2;++end)for(auto link:f.graph.nodes[i].ends[end])n.ends[end^flip[i]].push_back({ids[link.to],uint8_t(link.toEnd^flip[link.to])});
  }
  check(tg.validate().empty(),"transformed graph invariant");auto te=ts::collectReadThreadEvidence(tg,reads({p}),eligible);
  check(one(te,{uint64_t(ids[0]*2+flip[0]),uint64_t(ids[2]*2+flip[2]),uint64_t(ids[3]*2+flip[3])}),"node-order and unitig-strand invariance");
 }
 Fixture f;
 auto weakEnd=collect(f,reads({}, {f.left.substr(405)+f.repeat+f.right.substr(0,1)}));check(weakEnd.routes.empty()&&weakEnd.stats.readsNoFlanks==1,"one destination anchor insufficient");
 auto weakStart=collect(f,reads({}, {f.left.substr(499)+f.repeat+f.right.substr(0,95)}));check(weakStart.routes.empty()&&weakStart.stats.readsNoFlanks==1,"one source anchor insufficient");
 auto invalid=ts::collectReadThreadEvidence(f.graph,reads({f.pair(251)}),{});check(invalid.routes.empty()&&invalid.stats.readsExamined==0,"invalid eligibility mask refused");
 // Two exact interiors sharing endpoints must remain distinct evidence routes.
 Fixture interior(false,31);interior.graph.nodes[1].deleted=true;interior.graph.nodes[1].ends[1].clear();
 interior.graph.nodes[2].ends[0]={{0,1}};interior.graph.nodes[2].ends[1]={{3,0}};
 interior.graph.nodes[4].seq=interior.repeat;interior.graph.nodes[4].seq[70]=interior.repeat[70]=='A'?'C':'A';
 interior.graph.nodes[4].ends[0].clear();interior.graph.nodes[4].ends[1].clear();interior.eligible[4]=0;
 interior.graph.nodes[3].ends[0]={{2,1}};edge(interior.graph,0,4);edge(interior.graph,4,3);
 check(interior.graph.validate().empty(),"two-interior fixture invariant");
 std::string a=interior.left.substr(405)+interior.repeat+interior.right.substr(0,16),b=interior.left.substr(405)+interior.graph.nodes[4].seq+interior.right.substr(0,16);
 auto separate=collect(interior,reads({}, {a,b}));check(separate.routes.size()==2&&separate.routes[0].oriented.front()==separate.routes[1].oriented.front()&&separate.routes[0].oriented.back()==separate.routes[1].oriented.back()&&separate.routes[0].oriented!=separate.routes[1].oriented,"endpoint-equal interior-distinct routes not collapsed");
 // An even-k palindromic seed does not determine an oriented placement.
 ts::UnitigGraph palindrome;palindrome.setK(4);palindrome.nodes.resize(2);
 palindrome.nodes[0].seq="ATATA";palindrome.nodes[1].seq="ATACGGTAA";edge(palindrome,0,1);
 auto pe=ts::collectReadThreadEvidence(palindrome,reads({}, {"ATATACGGTAA"}),{1,1});
 check(pe.routes.empty(),"palindromic source seeds cannot establish strand");
 std::printf("%d checks passed\n",checks);
 return 0;
}
