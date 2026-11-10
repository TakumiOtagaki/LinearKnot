// Fixed CONTRAfold-v2 scores distributed by LinearSampling. All values are
// dimensionless log potentials. Pair terms are owned by the emission consumer.
#include "contra_scores.hh"
#include <cmath>
#define base_pair contrafold_base_pair
#include "utility.h"
#undef base_pair
namespace pkfit {
static int b(const std::string&s,int i){return i>=1 && i<=int(s.size())?GET_ACGU_NUM(s[i-1]):4;}
double pair(const std::string&s,int i,int j){return base_pair_score(b(s,i),b(s,j));}
double hairpin(const std::string&s,int i,int j){
  return score_hairpin(i-1,j-1,b(s,i),b(s,i+1),b(s,j-1),b(s,j));
}
double internal(const std::string&s,int i,int j,int k,int l){
  int u=k-i-1,v=j-l-1;
  if(u<0 || v<0 || u+v>30)return -INFINITY;
  if(u==0 && v==0)return helix_stacking_score(b(s,i),b(s,k),b(s,l),b(s,j));
  // Length tables have already accumulated their at-least features upstream.
  double value=0.;
  if(u==0 || v==0)value=bulge_length[u+v];
  else {
    value=internal_length[u+v];
    if(u<=4 && v<=4)value+=internal_explicit[u<=v?u*4+v:v*4+u];
    if(u==v)value+=internal_symmetric_length[std::min(u,15)];
    else value+=internal_asymmetry[std::min(std::abs(u-v),28)];
  }
  return value+score_junction_B(i-1,j-1,b(s,i),b(s,i+1),b(s,j-1),b(s,j))
    +score_junction_B(l-1,k-1,b(s,l),b(s,l+1),b(s,k-1),b(s,k))
    +score_single_nuc(i-1,j-1,k-1,l-1,b(s,k-1),b(s,l+1));
}
double stem(const std::string&s,int i,int j,bool multi,bool closing){
  double value=closing?
    score_junction_A(i-1,j-1,b(s,i),b(s,i+1),b(s,j-1),b(s,j),s.size()):
    score_junction_A(j-1,i-1,b(s,j),b(s,j+1),b(s,i-1),b(s,i),s.size());
  return value+(multi?multi_paired:external_paired);
}
double external_unpaired_score(){return external_unpaired;}
double multi_unpaired_score(){return multi_unpaired;}
double multi_base_score(){return multi_base;}
double multi_paired_score(){return multi_paired;}
}
