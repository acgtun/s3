#ifndef EVALUATION_H_
#define EVALUATION_H_

#include "bio_util.h"
#include "option.h"

#include <vector>
#include <fstream>
#include <tr1/unordered_map>

using std::cout;
using std::endl;
using std::vector;
using std::tr1::unordered_map;

typedef unordered_map<string, vector<M8Results> > ResultsType;

class Evaluation {
 public:
  Evaluation(const string& _results_file, Evaluation* _ground_truth)
      : results_file(_results_file),
        ground_truth(_ground_truth) {
    Option::GetOption("-evalue", evalue_threshold, 10);
    Option::ChkStrExist("-sort", is_sort);
    Option::GetOption("-K", K, 1);
    Option::GetOption("-M", M, 10);
  }
  ~Evaluation() {
  }

  void ReadResults();
  void ShowReadResults(ofstream& fout);
  void CompareTopKHitProtein(ofstream& fout);
  void CompareTopKHitProtein_INFO();

  string results_file;
  ResultsType search_results;

  double evalue_threshold;
  bool is_sort;
  uint32_t K;
  uint32_t M;
  Evaluation* ground_truth;
};

#endif /* EVALUATION_H_ */
