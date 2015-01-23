#ifndef LOCALALIGNMENT_H_
#define LOCALALIGNMENT_H_

#include "bio_util.h"
#include "evalue.h"

#include <vector>
using std::vector;

namespace local_alignment {

#define DIAG ('a')
#define UP ('b')
#define LEFT ('c')
#define STOPorSTART ('d')

class LocalAlignment {
 public:
  LocalAlignment(const Evalue* _evalue, const uint32_t& _max_rows,
                 const uint32_t& _max_cols);
  ~LocalAlignment();

  bool RunLocalAlignment(const string& U, const vector<char>& V,
                         M8Results& res);
 private:
  vector<char> rU;
  vector<char> rV;
  vector<char> midline;
  //char * rU;
  //char * rV;
  //char * midline;
  int rs;

  clock_t start_t;
  uint64_t sum_time;

  vector<vector<int> > s;
  vector<vector<char> > l;

  vector<vector<int> > g;
  vector<vector<int> > h;

  uint32_t n;
  uint32_t m;

  int alignScore;
  int nIdentity;
  double e_value;

  int gapopen;
  int gapextension;

  const Evalue* evalue;

  uint32_t max_rows;
  uint32_t max_cols;

  void DisplayAlignment();
  void stringReverse(vector<char>& str, const int & n);
  int MaxOfFour(const int & s1, const int & s2, const int & s3, const int & s4);
  char Direction(const int & s1, const int & s2, const int & s3,
                 const int & s4);
};

}  // namespace local_alignment

#endif /* LOCALALIGNMENT_H_ */
