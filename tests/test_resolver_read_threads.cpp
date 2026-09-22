// Actual anchoring/support/reciprocal-resolution regression for exact read paths.
#define main helper_unit_main_not_called
#include "test_read_thread_evidence.cpp"
#undef main
#include "resolve.h"
#include <cstring>
namespace {
struct Resolved {
 std::vector<std::string> contigs;std::vector<ts::ResolvedPath> paths;
 ts::ResolveStats stats;std::string log;
};
Resolved resolve(const ts::UnitigGraph& graph,const ts::SequenceStore& input,const char* flag){
 if(flag)setenv("TESSERACT_EXACT_READ_THREADS",flag,1);else unsetenv("TESSERACT_EXACT_READ_THREADS");
 FILE* capture=tmpfile();check(capture!=nullptr,"stderr capture");int saved=dup(2);check(saved>=0,"stderr descriptor");
 fflush(stderr);dup2(fileno(capture),2);
 ts::PairedResolver resolver(graph,input,1,2,1.02,0.0);resolver.buildSupport();
 Resolved out;std::vector<double> coverage;resolver.resolve(out.contigs,coverage);out.paths=resolver.paths();out.stats=resolver.stats();
 fflush(stderr);dup2(saved,2);close(saved);rewind(capture);char buffer[4096];size_t n=0;
 while((n=fread(buffer,1,sizeof(buffer),capture))) { out.log.append(buffer,n); }
 fclose(capture);return out;
}
bool hasPath(const Resolved& r,ts::ReadThreadPath path){
 for(const auto& p:r.paths) { if(p.oriented==path||p.oriented==ts::reverseReadThreadPath(path))return true; }
 return false;
}
bool includesPath(const Resolved& r,const ts::ReadThreadPath& path){
 for(const auto& p:r.paths){auto rev=ts::reverseReadThreadPath(p.oriented);
  if(std::search(p.oriented.begin(),p.oriented.end(),path.begin(),path.end())!=p.oriented.end()||
     std::search(rev.begin(),rev.end(),path.begin(),path.end())!=rev.end())return true;
 }return false;
}
Fixture interiorFixture(){
 Fixture f(false,31);f.graph.nodes[1].deleted=true;f.graph.nodes[1].ends[1].clear();
 f.graph.nodes[2].ends[0]={{0,1}};f.graph.nodes[2].ends[1]={{3,0}};f.graph.nodes[2].coverage=50;
 f.graph.nodes[4].seq=f.repeat;f.graph.nodes[4].seq[70]=f.repeat[70]=='A'?'C':'A';f.graph.nodes[4].coverage=40;
 f.graph.nodes[4].ends[0].clear();f.graph.nodes[4].ends[1].clear();f.eligible[4]=0;
 f.graph.nodes[3].ends[0]={{2,1}};edge(f.graph,0,4);edge(f.graph,4,3);
 ts::Unitig other;other.seq=f.repeat.substr(110)+dna(500,155);other.coverage=20;f.graph.nodes.push_back(other);f.eligible.push_back(1);
 edge(f.graph,2,5);edge(f.graph,4,5);check(f.graph.validate().empty(),"interior ambiguity graph");return f;
}
}
int main(int argc,char** argv){
 setenv("TESSERACT_COMMON_PREFIX","0",1);setenv("TESSERACT_JOIN_TRACE","1",1);
 for(auto flag:{"TESSERACT_ROUTE_DISTANCE","TESSERACT_WEIGHTED_RESOLVER_COVERAGE","TESSERACT_WEIGHTED_ELIGIBLE_COVERAGE","TESSERACT_PREFIX_SNP_BUBBLES","TESSERACT_EXCLUDE_SHARED_REPEAT_SUPPORT","TESSERACT_SHARED_SUPPORT_AUDIT","TESSERACT_REQUIRE_SUPPORT_SINGLE"})unsetenv(flag);
 if(argc==2&&std::string(argv[1])=="--dump-control"){
  for(int length:{251,301}){
   Fixture f;auto pair=f.pair(length),alternative=f.pair(length,true);
   std::vector<std::vector<std::pair<std::string,std::string>>> cases{{pair},{pair,pair},{pair,pair,pair,alternative}};
   auto supported=cases[1];for(int i=0;i<5;++i)supported.push_back({f.graph.nodes[0].seq.substr(100+i,70),ts::reverseComplement(f.graph.nodes[4].seq.substr(200+i,70))});cases.push_back(supported);
   for(size_t i=0;i<cases.size();++i){auto result=resolve(f.graph,reads(cases[i]),nullptr);
    std::printf("CASE %d %zu links=%zu joined=%zu\n",length,i,result.stats.pairsLinking,result.stats.unitigsJoined);
    for(const auto& contig:result.contigs) { std::printf("%s\n",contig.c_str()); }
    std::printf("%s",result.log.c_str());
   }
  }
  return 0;
 }
 for(int length:{251,301}){
  Fixture f;auto pair=f.pair(length);auto input=reads({pair,pair});
  auto baseline=resolve(f.graph,input,nullptr),off=resolve(f.graph,input,"0"),nonnumeric=resolve(f.graph,input,"true"),candidate=resolve(f.graph,input,"1");
  // 1.3.0: exact threads are ON unless TESSERACT_EXACT_READ_THREADS=0.
  check(baseline.contigs==candidate.contigs&&baseline.contigs==nonnumeric.contigs,"default-on: unset and nonzero flag identity");
  check(baseline.log==candidate.log&&baseline.log==nonnumeric.log,"default-on trace identity");
  check(off.stats.pairsLinking==0,"fixture loses all legacy inter-unitig pair support");
  auto reversedGraph=f.graph;for(auto& node:reversedGraph.nodes){node.seq=ts::reverseComplement(node.seq);std::swap(node.ends[0],node.ends[1]);for(auto& links:node.ends)for(auto& link:links)link.toEnd^=1;}
  check(reversedGraph.validate().empty(),"reverse represented resolver graph");
  check(hasPath(resolve(reversedGraph,input,"1"),{1,5,7}),"full resolver preserves reverse-unitig orientation");
  check(!hasPath(off,{0,4,6}),"=0 legacy resolver remains unresolved");
  check(hasPath(candidate,{0,4,6}),"exact-thread fallback creates correct observed join");
  check(candidate.log.find("molecules=2 fresh=2 pairedExcluded=0")!=std::string::npos,"two physical molecules, four reads, two fresh votes");
  const std::string spelled=f.graph.nodes[0].seq+f.graph.nodes[2].seq.substr(126)+f.graph.nodes[3].seq.substr(126);
  bool correctSequence=false;for(const auto& s:candidate.contigs)correctSequence|=s==spelled||ts::reverseComplement(s)==spelled;
  check(correctSequence,"joined FASTA exactly spells observed graph connector");
  auto oneMolecule=resolve(f.graph,reads({pair}),"1");check(!hasPath(oneMolecule,{0,4,6}),"one molecule with two matching mates cannot meet two-vote floor");
  auto alternative=f.pair(length,true);auto conflict=resolve(f.graph,reads({pair,pair,pair,alternative}),"1");
  check(!hasPath(conflict,{0,4,6})&&!hasPath(conflict,{0,4,8}),"minority contrary endpoint vetoes fallback");
  auto mateConflict=resolve(f.graph,reads({{pair.first,alternative.second},{pair.first,alternative.second}}),"1");
  check(!hasPath(mateConflict,{0,4,6})&&!hasPath(mateConflict,{0,4,8})&&mateConflict.log.find("conflictingMates=2")!=std::string::npos,"full resolver refuses inconsistent mates");
  Fixture amb(true);auto ambiguous=resolve(amb.graph,input,"1");check(!hasPath(ambiguous,{0,4,6})&&!hasPath(ambiguous,{0,4,8}),"ambiguous complete read mappings refuse fallback");
  auto masked=reads({pair,pair});for(size_t r=0;r<masked.size();++r)masked.maskRange(r,100,101);
  auto unknown=resolve(f.graph,masked,"1");check(!hasPath(unknown,{0,4,6}),"corrector masks cannot be treated as raw sequence evidence");
  auto missingMate=resolve(f.graph,reads({{pair.first,std::string(size_t(length),'N')},{pair.first,std::string(size_t(length),'N')}}),"1");
  check(hasPath(missingMate,{0,4,6})&&missingMate.stats.pairsLinking==0,"two valid reads with unavailable mates can phase observed connector");
  std::vector<std::pair<std::string,std::string>> supported{pair,pair};
  for(int i=0;i<5;++i)supported.push_back({f.graph.nodes[0].seq.substr(100+i,70),ts::reverseComplement(f.graph.nodes[4].seq.substr(200+i,70))});
  auto chosen=resolve(f.graph,reads(supported),"0"),kept=resolve(f.graph,reads(supported),"1");
  check(hasPath(chosen,{0,4,8}),"legacy fixture has supported competing choice");
  check(chosen.contigs==kept.contigs,"existing accepted paired choice unchanged");
  // Fresh exact evidence is not a second use of already linked endpoint pairs.
  std::string fragment=f.left.substr(399)+f.repeat+f.right.substr(0,80);
  const std::pair<std::string,std::string> linked{fragment.substr(0,251),ts::reverseComplement(fragment.substr(70,251))};
  auto reused=resolve(f.graph,reads({linked}),"1");
  check(reused.stats.pairsLinking==1,"exclusion fixture has existing paired endpoint");
  check(reused.log.find("fresh=0 pairedExcluded=1")!=std::string::npos,"already paired molecule excluded from exact route count");
 }
 // Keep initial production scope paired; standalone helper single-end tests remain separate.
 Fixture f;auto p=f.pair(301);auto single=reads({}, {p.first,p.first});
 check(resolve(f.graph,single,"0").contigs==resolve(f.graph,single,"1").contigs,"pure single-end resolver behavior unchanged");
 // A new join from A's right port must not consume B before B's supported
 // right-port join to Z. The trace verifies legacy arbitration runs first.
 Fixture priority;ts::Unitig z;z.seq=priority.graph.nodes[3].seq.substr(500)+dna(500,91);z.coverage=20;
 priority.graph.nodes.push_back(z);edge(priority.graph,3,5);
 auto pp=priority.pair(301);std::vector<std::pair<std::string,std::string>> molecules{pp,pp};
 for(int i=0;i<5;++i)molecules.push_back({priority.graph.nodes[3].seq.substr(200+i,70),ts::reverseComplement(priority.graph.nodes[5].seq.substr(200+i,70))});
 auto pr=resolve(priority.graph,reads(molecules),"1");size_t first=pr.log.find("[acceptedjoin]");
 check(first!=std::string::npos,"priority fixture joins");std::string event=pr.log.substr(first,pr.log.find('\n',first)-first);
 check(event.find("route=6,10")!=std::string::npos,"existing opposite-port join is accepted before fallback consumes chain");
 check(includesPath(pr,{6,10})&&includesPath(pr,{0,4,6}),"existing join and new exact join both retained");
 // Complete connector reciprocity, not mere endpoint agreement.
 auto interior=interiorFixture();std::string r1=interior.left.substr(405)+interior.repeat+interior.right.substr(0,16);
 std::string r2=interior.left.substr(405)+interior.graph.nodes[4].seq+interior.right.substr(0,16);
 auto pair1=std::make_pair(r1,ts::reverseComplement(r1)),pair2=std::make_pair(r2,ts::reverseComplement(r2));
 auto competingRoutes=resolve(interior.graph,reads({pair1,pair1,pair2,pair2}),"1");
 check(!includesPath(competingRoutes,{0,4,6})&&!includesPath(competingRoutes,{0,8,6}),"distinct interiors sharing endpoints veto exact fallback");
 auto mismatch=resolve(interior.graph,reads({pair2,pair2}),"1");
 check(mismatch.log.find("nominations=1")!=std::string::npos,"mismatch fixture nominates only unsupported direction");
 check(!includesPath(mismatch,{0,8,6})&&!includesPath(mismatch,{0,4,6}),"opposite legacy connector disagreement refuses exact fallback join");
 std::printf("%d full resolver exact-thread checks passed\n",checks);
 return 0;
}
