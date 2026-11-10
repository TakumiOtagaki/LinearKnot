// Independent tree scorer for PK-free oracle checks, without the DP grammar.
#include "contra_scores.hh"
#include <iostream>
#include <iomanip>
#include <vector>
#include <functional>
int main(){std::string s,g;std::cout<<std::setprecision(17);
 while(std::cin>>s>>g){
  if(s.size()!=g.size())return 2;
  std::vector<int> p(s.size()+1,-1),stack;
  for(int i=1;i<=int(s.size());++i){if(g[i-1]=='(')stack.push_back(i);
    else if(g[i-1]==')'){if(stack.empty())return 2;int j=stack.back();stack.pop_back();p[i]=j;p[j]=i;}}
  if(!stack.empty())return 2;
  std::function<double(int,int)> loop=[&](int i,int j){
    std::vector<std::pair<int,int>> children;int unpaired=0;
    for(int k=i+1;k<j;){if(p[k]>k){children.emplace_back(k,p[k]);k=p[k]+1;}else{++unpaired;++k;}}
    double score=pkfit::pair(s,i,j);
    if(children.empty())return score+pkfit::hairpin(s,i,j);
    if(children.size()==1){auto[k,l]=children[0];return score+pkfit::internal(s,i,j,k,l)+loop(k,l);}
    score+=pkfit::multi_base_score()+pkfit::stem(s,i,j,true,true)+unpaired*pkfit::multi_unpaired_score();
    for(auto[k,l]:children)score+=loop(k,l)+pkfit::stem(s,k,l,true,false);
    return score;
  };
  double score=0.;for(int i=1;i<=int(s.size());){if(p[i]>i){score+=loop(i,p[i])+pkfit::stem(s,i,p[i],false,false);i=p[i]+1;}
    else {score+=pkfit::external_unpaired_score();++i;}}
  std::cout<<score<<'\n';
 }
}
