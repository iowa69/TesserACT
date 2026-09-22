#include "graph_coverage.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <random>

namespace {
int checks=0;
void check(bool ok,const char* message){++checks;if(!ok){std::fprintf(stderr,"FAIL: %s\n",message);std::exit(1);}}
void add(ts::UnitigGraph& graph,size_t mass,double coverage,bool deleted=false){
    ts::Unitig u;u.seq.assign(mass+graph.k()-1,'A');u.coverage=coverage;u.deleted=deleted;graph.nodes.push_back(std::move(u));
}
}
int main(){
    ts::UnitigGraph graph;graph.setK(127);
    check(ts::graphLengthWeightedMedianCoverage(graph)==0,"empty graph");
    add(graph,4000,30);add(graph,500,90);add(graph,500,4);
    check(ts::graphLengthWeightedMedianCoverage(graph)==30,"minority repeat and error mass do not set median");
    add(graph,1000000,0);add(graph,1000000,1,true);
    add(graph,1000000,std::numeric_limits<double>::infinity());
    add(graph,1000000,std::numeric_limits<double>::quiet_NaN());
    check(ts::graphLengthWeightedMedianCoverage(graph)==30,"zero nonfinite and deleted nodes carry no weight");
    ts::Unitig tiny;tiny.seq.assign(126,'A');tiny.coverage=1;graph.nodes.push_back(tiny);
    check(ts::graphLengthWeightedMedianCoverage(graph)==30,"shorter-than-k nodes carry no weight");
    std::mt19937 random(56191);
    for(int attempt=0;attempt<100;++attempt){
        ts::UnitigGraph split;split.setK(127);
        for(auto original: std::vector<std::pair<size_t,double>>{{4000,30},{500,90},{500,4}}){
            size_t remaining=original.first;
            while(remaining){size_t part=std::min(remaining,size_t(1+random()%100));add(split,part,original.second);remaining-=part;}
        }
        // Every piece is <2k: the estimator must still represent all 5000
        // independent graph bases. Full stored-sequence length is not mass.
        check(ts::graphLengthWeightedMedianCoverage(split)==30,"arbitrary unitig splitting preserves estimate");
        std::shuffle(split.nodes.begin(),split.nodes.end(),random);
        check(ts::graphLengthWeightedMedianCoverage(split)==30,"node-order invariant");
    }
    ts::UnitigGraph fragmented;fragmented.setK(127);
    add(fragmented,10000,30);
    for(int i=0;i<200;++i)add(fragmented,10,3);
    check(ts::graphLengthWeightedMedianCoverage(fragmented)==30,"many tiny low-depth fragments cannot outvote chromosome mass");
    ts::UnitigGraph overlap;overlap.setK(127);
    add(overlap,1100,30);
    for(int i=0;i<100;++i)add(overlap,10,3);
    check(ts::graphLengthWeightedMedianCoverage(overlap)==30,"duplicated k-1 boundaries do not inflate fragment weight");
    ts::UnitigGraph empty;empty.setK(127);add(empty,100,0);add(empty,100,30,true);
    check(ts::graphLengthWeightedMedianCoverage(empty)==0,"no positive live mass");
    ts::UnitigGraph equal;equal.setK(31);add(equal,100,10);add(equal,100,20);
    check(ts::graphLengthWeightedMedianCoverage(equal)==10,"exact half-mass tie selects lower weighted median");
    // Conditional estimator preserves the original >=2k population.
    check(ts::graphEligibleLengthWeightedMedianCoverage(graph)==30,"eligible estimator ignores invalid/deleted depths");
    check(ts::graphEligibleLengthWeightedMedianCoverage(empty)==0,"eligible empty population returns legacy-fallback sentinel");
    ts::UnitigGraph edge;edge.setK(31);
    add(edge,31,2);  // len = 2k-1, excluded
    check(ts::graphEligibleLengthWeightedMedianCoverage(edge)==0,"below 2k eligibility boundary");
    add(edge,32,18); // len = 2k, included
    check(ts::graphEligibleLengthWeightedMedianCoverage(edge)==18,"exactly 2k eligibility boundary");
    ts::UnitigGraph stress;stress.setK(99);
    add(stress,3000,18);
    for(int i=0;i<500;++i)add(stress,10,2);
    check(ts::graphLengthWeightedMedianCoverage(stress)==2,"all-node experiment preserves its known short-mass behavior");
    check(ts::graphEligibleLengthWeightedMedianCoverage(stress)==18,"eligible population excludes dominating short low-depth mass");
    ts::UnitigGraph biased;biased.setK(31);add(biased,5000,30);add(biased,500,90);
    for(int i=0;i<50;++i)add(biased,32,4);
    check(ts::graphEligibleLengthWeightedMedianCoverage(biased)==30,"many eligible short pieces cannot outvote chromosome mass");
    for(int attempt=0;attempt<100;++attempt){
        ts::UnitigGraph split;split.setK(127);
        for(auto original: std::vector<std::pair<size_t,double>>{{4000,30},{500,90},{500,4}}){
            size_t remaining=original.first;
            while(remaining>=256){
                const size_t largest=std::min(size_t(400),remaining-128);
                const size_t part=128+random()%(largest-128+1);
                add(split,part,original.second);remaining-=part;
            }
            add(split,remaining,original.second);
        }
        check(ts::graphEligibleLengthWeightedMedianCoverage(split)==30,"eligible-preserving splits leave conditional median unchanged");
        std::shuffle(split.nodes.begin(),split.nodes.end(),random);
        check(ts::graphEligibleLengthWeightedMedianCoverage(split)==30,"conditional estimator is node-order invariant");
    }
    check(ts::graphEligibleLengthWeightedMedianCoverage(equal)==10,"eligible weighted half-mass tie uses lower median");
    // ANTI-PATTERN GUARD (directive S25): this estimator must be a weighted
    // MEDIAN, never a weighted MEAN. A single very deep node must not drag the
    // estimate the way an average would.
    {
        ts::UnitigGraph mm;mm.setK(127);
        add(mm,1000,10);add(mm,1000,10);add(mm,1000,1000);
        const double got=ts::graphLengthWeightedMedianCoverage(mm);
        // weighted mean would be (10+10+1000)/3 = 340 at equal mass
        check(got==10,"weighted median ignores a lone deep node");
        check(std::fabs(got-340.0)>1.0,"estimator is a median, NOT a weighted mean");
    }
    // Reproduces the measured pathology on senterica/GCA006365335v1: a cloud of
    // short unitigs sitting at the carry floor outvotes the chromosome by COUNT
    // while carrying almost none of it by LENGTH. Unweighted node median = 8;
    // length-weighted median must recover the true ~35x depth.
    {
        ts::UnitigGraph path;path.setK(99);
        for(int i=0;i<327;++i) add(path,50,4.0);      // carry-floor cloud, 2% of mass
        for(int i=0;i<435;++i) add(path,60,8.0);      // low-depth fragments, 3% of mass
        for(int i=0;i<192;++i) add(path,20000,35.0);  // the actual chromosome, 94%
        const double got=ts::graphLengthWeightedMedianCoverage(path);
        check(got==35.0,"carry-floor cloud does not set the baseline (senterica shape)");
        check(got>1.6*8.0,"recovered depth clears the repeat threshold the old median produced");
    }
    std::printf("graph coverage: %d checks passed\n",checks);
}
