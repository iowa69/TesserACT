// Directly check the constructor's experiment selection without exposing a
// production accessor solely for tests. Include dependencies before the macro.
#include "graph.h"
#include "seqio.h"
#include <unordered_map>
#define private public
#include "resolve.h"
#undef private
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>

namespace {
int checks=0;
void check(bool ok,const char* message){++checks;if(!ok){std::fprintf(stderr,"FAIL: %s\n",message);std::exit(1);}}
void add(ts::UnitigGraph& graph,size_t mass,double coverage){
    ts::Unitig u;u.seq.assign(mass+graph.k()-1,'A');u.coverage=coverage;graph.nodes.push_back(std::move(u));
}
double selected(const ts::UnitigGraph& graph,const char* all,const char* eligible){
    if(all)setenv("TESSERACT_WEIGHTED_RESOLVER_COVERAGE",all,1);else unsetenv("TESSERACT_WEIGHTED_RESOLVER_COVERAGE");
    if(eligible)setenv("TESSERACT_WEIGHTED_ELIGIBLE_COVERAGE",eligible,1);else unsetenv("TESSERACT_WEIGHTED_ELIGIBLE_COVERAGE");
    ts::SequenceStore reads;
    ts::PairedResolver resolver(graph,reads,1,3,1.5);
    return resolver.medianCoverage_;
}
}
int main(){
    ts::UnitigGraph stress;stress.setK(99);add(stress,3000,18);
    for(int i=0;i<500;++i)add(stress,10,2);
    check(selected(stress,nullptr,nullptr)==18,"no flag preserves legacy");
    check(selected(stress,"1",nullptr)==2,"all-node experiment retained");
    check(selected(stress,nullptr,"1")==18,"eligible experiment excludes short mass");
    check(selected(stress,"1","1")==18,"eligible experiment has explicit precedence");
    check(selected(stress,"1","0")==2,"zero eligible flag leaves all-node selection");
    check(selected(stress,"1","true")==2,"only literal 1 enables eligible flag");
    check(selected(stress,"true",nullptr)==18,"only literal 1 enables all-node flag");
    ts::UnitigGraph empty;empty.setK(99);add(empty,10,2);
    check(selected(empty,nullptr,"1")==0,"no eligible population retains legacy zero");
    check(selected(empty,"1","1")==0,"eligible empty fallback does not switch to all-node estimator");
    ts::UnitigGraph negative;negative.setK(99);add(negative,3000,-4);
    check(selected(negative,nullptr,"1")==-4,"no valid mass retains legacy negative value exactly");
    ts::UnitigGraph nan;nan.setK(99);add(nan,3000,std::numeric_limits<double>::quiet_NaN());
    check(std::isnan(selected(nan,nullptr,"1")),"no valid mass retains legacy NaN");
    ts::UnitigGraph biased;biased.setK(31);add(biased,5000,30);add(biased,500,90);
    for(int i=0;i<50;++i)add(biased,32,4);
    check(selected(biased,nullptr,nullptr)==4,"legacy one-vote-per-unitig behavior preserved");
    check(selected(biased,nullptr,"1")==30,"eligible mode changes weights within original length population");
    unsetenv("TESSERACT_WEIGHTED_RESOLVER_COVERAGE");unsetenv("TESSERACT_WEIGHTED_ELIGIBLE_COVERAGE");
    std::printf("resolver coverage modes: %d checks passed\n",checks);
}
