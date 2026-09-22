#include "polish_quality.h"
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <random>
#include <vector>
using namespace ts::quality_consensus;
namespace {
size_t checks = 0;
void check(bool ok, const char* label) { ++checks; if (!ok) { std::fprintf(stderr,"FAIL: %s\n",label); std::exit(1); } }
void add(Pileup& p, int base, int q, int strand) { Observation o{uint8_t(base),uint8_t(q),uint8_t(strand)}; addFragment(p,&o,1); }
}
int main() {
    check(logOdds(0)==0 && logOdds(1)==0,"Q0 and below-chance Q1 are uniform");
    for(unsigned q=2;q<=40;++q) check(logOdds(q)>logOdds(q-1),"informative odds increase with Q");
    for(unsigned q=40;q<=93;++q) check(logOdds(q)==logOdds(40),"fixed Q40 cap");
    for(int b=0;b<4;++b) {
        Pileup single,pair;
        for(int i=0;i<15;++i) { add(single,b,30,i%2); Observation os[2]={{uint8_t(b),30,uint8_t(i%2)},{uint8_t(b),30,uint8_t(i%2)}};addFragment(pair,os,2); }
        check(single.scores==pair.scores && single.fragments==pair.fragments,"identical overlapping mates do not double evidence/depth");
        check(decide(single,(b+1)%4,15,.9).reason==DecisionReason::Accept,"15 fragments both orientations accepted");
        check(decide(single,b,15,.9).reason==DecisionReason::CurrentAllele,"current allele preserved");
        check(decide(single,(b+1)%4,16,.9).reason==DecisionReason::BelowDepth,"existing depth floor inclusive");
    }
    Pileup neutral; Observation zs[2]={{0,0,0},{1,0,1}};
    for(int i=0;i<100;++i) check(addFragment(neutral,zs,2)==AddResult::Uniform,"Q0 pair uniform");
    check(neutral.fragments==0 && neutral.orientations[0]==0 && neutral.observations==200,"Q0 adds neither informative depth nor strand support");
    Pileup half;Observation hs[2]={{0,40,0},{1,0,1}};addFragment(half,hs,2);
    check(half.scores[0]==logOdds(40) && half.scores[1]==0 && half.fragments==1,"Q0 remains in pair mean denominator");
    check(half.orientations[0]==1 && half.orientations[1]==0,"uniform mate cannot supply strand support");
    Pileup conflict;Observation cs[2]={{0,35,0},{1,35,1}};
    for(int i=0;i<20;++i)addFragment(conflict,cs,2);
    check(conflict.fragments==20 && conflict.scores[0]==conflict.scores[1],"conflicting mates average and depth once");
    check(decide(conflict,2,15,.1).reason==DecisionReason::Tie,"no alphabet-order tie choice even at permissive fraction");
    Pileup low;
    for(int i=0;i<14;++i)add(low,0,40,i%2);
    for(int i=0;i<100;++i)addFragment(low,zs,2);
    check(decide(low,1,15,.9).reason==DecisionReason::BelowDepth,"14 informative plus many Q0 stays below15");
    Pileup one;
    for(int i=0;i<20;++i)add(one,0,30,0);
    check(decide(one,1,15,.9).reason==DecisionReason::OneOrientation,"one strand abstains despite high score");
    // Both orientations may be supplied by a single pair; this is a safeguard,
    // not a requirement for independent fragments on each orientation.
    Observation both[2]={{0,30,0},{0,30,1}}; addFragment(one,both,2);
    check(decide(one,1,15,.9).reason==DecisionReason::Accept,"same fragment may supply both orientation bits");
    Pileup weak;
    for(int i=0;i<8;++i)add(weak,0,2,i%2);
    for(int i=0;i<7;++i)add(weak,1,2,i%2);
    check(decide(weak,1,15,.9).reason==DecisionReason::BelowPosterior,"adequate depth and both strands cannot replace confidence floor");
    auto boundary=decide(weak,1,15,0);
    check(decide(weak,1,15,boundary.posterior).reason==DecisionReason::Accept,"score-derived fraction threshold inclusive");
    check(decide(weak,1,15,std::nextafter(boundary.posterior,1.0)).reason==DecisionReason::BelowPosterior,"above exact score threshold rejected");
    check(decide(weak,1,15,std::numeric_limits<double>::quiet_NaN()).reason==DecisionReason::InvalidThreshold,"nonfinite threshold rejected");
    Pileup overflow;overflow.scores[0]=std::numeric_limits<uint64_t>::max();auto before=overflow.scores;
    Observation good{0,40,0};check(addFragment(overflow,&good,1)==AddResult::Overflow && overflow.scores==before,"score overflow atomic fail closed");
    check(decide(overflow,1,15,.9).reason==DecisionReason::Overflow,"overflow cannot propose");
    overflow={};overflow.fragments=std::numeric_limits<uint64_t>::max();check(addFragment(overflow,&good,1)==AddResult::Overflow,"fragment count overflow");
    overflow={};overflow.observations=std::numeric_limits<uint64_t>::max();check(addFragment(overflow,&good,1)==AddResult::Overflow,"observation count overflow");
    Pileup invalid;Observation bad{4,20,0};check(addFragment(invalid,&bad,1)==AddResult::Invalid && invalid.fragments==0,"invalid allele rejected");
    bad={0,255,0};check(addFragment(invalid,&bad,1)==AddResult::Invalid,"missing quality never capped into strong evidence");
    check(addFragment(invalid,&good,3)==AddResult::Invalid,"more than two observations requires caller deduplication");
    // Integer addition is exactly invariant to fragment/mate order and allele
    // relabeling. Compare posterior to an independent unrounded log model.
    std::mt19937 rng(92026);
    for(int round=0;round<500;++round) {
        std::vector<std::vector<Observation>> fragments;
        Pileup p,reverse,renamed;double reference[4]={};
        for(int f=0;f<40;++f) {
            std::vector<Observation> os;
            for(unsigned j=0,n=1+rng()%2;j<n;++j)os.push_back({uint8_t(rng()%4),uint8_t(rng()%94),uint8_t(rng()%2)});
            fragments.push_back(os);addFragment(p,os.data(),os.size());
            for(auto o:os) {
                const double e=std::min(.75,std::pow(10.0,-double(std::min<unsigned>(o.phred,40))/10));
                reference[o.base]+=std::log(3*(1-e)/e)/os.size();
                o.base=3-o.base;o.orientation=1-o.orientation;
            }
        }
        std::shuffle(fragments.begin(),fragments.end(),rng);
        for(auto os:fragments) {
            std::reverse(os.begin(),os.end());addFragment(reverse,os.data(),os.size());
            for(auto& o:os){o.base=3-o.base;o.orientation=1-o.orientation;}addFragment(renamed,os.data(),os.size());
        }
        check(p.scores==reverse.scores && p.fragments==reverse.fragments && p.orientations==reverse.orientations,"fragment/mate order invariance");
        for(int b=0;b<4;++b)check(p.scores[b]==renamed.scores[3-b],"reverse-complement allele symmetry");
        double maximum=*std::max_element(reference,reference+4),denom=0;
        for(double v:reference)denom+=std::exp(v-maximum);
        check(std::abs(decide(p,-1,15,.9).posterior-1/denom)<0.00001,"fixed point matches independent Phred likelihood within rounding error");
    }
    std::printf("%zu quality helper checks passed\n",checks);
}
