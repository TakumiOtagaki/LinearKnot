#pragma once
#include <string>
namespace pkfit {
inline int pair_kind(char a,char b){
  if((a=='A'&&b=='U')||(a=='U'&&b=='A'))return 0;
  if((a=='C'&&b=='G')||(a=='G'&&b=='C'))return 1;
  if((a=='G'&&b=='U')||(a=='U'&&b=='G'))return 2;
  return -1;
}
double hairpin(const std::string&, int, int);
double internal(const std::string&, int, int, int, int);
double stem(const std::string&, int, int, bool, bool);
double pair(const std::string&, int, int);
double external_unpaired_score();
double multi_unpaired_score();
double multi_base_score();
double multi_paired_score();
}
