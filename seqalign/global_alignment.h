#ifndef GLOCALALIGNMENT_H_
#define GLOCALALIGNMENT_H_

#include "./../util/bio_util.h"

namespace global_alignment {

class GlobalAlignment {
 public:
  GlobalAlignment();
  ~GlobalAlignment();
  int RunGlobalAlignment(const char * U, const char * V, const usint32_t & n,
                         const usint32_t & m);

 private:
//  clock_t start_t_;
//  usint64_t sum_time_;

  vector<vector<int> > s;
  vector<vector<int> > g;
  vector<vector<int> > h;

  int gapopen_;
  int gapextension_;

};

}  // namespace global_alignment

#endif /* GLOCALALIGNMENT_H_ */
