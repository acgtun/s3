#ifndef LOCALALIGNMENT_H_
#define LOCALALIGNMENT_H_

#include "./../util/bio_util.h"
#include "./../util/evalue.h"
#include <time.h>

namespace local_alignment {

//#define I (i - 1)
//#define J (j - 1)
//#define P (p - 1)
//#define Q (q - 1)

#define DIAG ('a')
#define UP ('b')
#define LEFT ('c')
#define STOPorSTART ('d')

struct AlignResult {
  string align_query;
  string midline;
  string align_sequence;
  int align_score;
  double e_value;
  double s_prime;
};

class LocalAlignment {
 public:
  LocalAlignment(const Evalue* evalue, const usint32_t& max_rows,
                 const usint32_t& max_cols);
  ~LocalAlignment();

  bool RunLocalAlignment(const char * U, const char * V, const usint32_t & Ul,
                         const usint32_t & Vl, M8Results& res);
 private:
  char * rU;
  char * rV;
  char * midline;
  int rs;

  clock_t start_t_;
  usint64_t sum_time_;

  vector<vector<int> > s;
  vector<vector<char> > l;

  vector<vector<int> > g;
  vector<vector<int> > h;

  usint32_t n;
  usint32_t m;

  int alignScore;
  int nIdentity;
  double e_value;

  int gapopen_;
  int gapextension_;

  const Evalue* evalue_;

  usint32_t max_rows_;
  usint32_t max_cols_;

  void DisplayAlignment();
  void stringReverse(char * str, const int & n);
  int MaxOfFour(const int & s1, const int & s2, const int & s3, const int & s4);
  char Direction(const int & s1, const int & s2, const int & s3,
                 const int & s4);
};

}  // namespace local_alignment

#endif /* LOCALALIGNMENT_H_ */
