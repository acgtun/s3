#ifndef EVALUATION_H_
#define EVALUATION_H_

#include "./../util/bio_util.h"
#include "./../util/option.h"

#include <map>

typedef map<string, vector<M8Results> > ResultsType;

class Evaluation {
 public:
  Evaluation(const string& results_file, Evaluation* ground_truth)
      : results_file_(results_file),
        ground_truth_(ground_truth) {
    Option::GetOption("-evalue", evalue_threshold_, 10);
    Option::ChkStrExist("-sort", is_sort_);
    Option::GetOption("-K", K_, 1);
    Option::GetOption("-M", M_, 10);
  }
  ~Evaluation() {

  }

  void ReadResults();
  void ShowReadResults(ofstream& fout);
  void CompareTopKHitProtein(ofstream& fout);
  void CompareTopKHitProtein_INFO();

  string results_file_;
  ResultsType search_results_;

  double evalue_threshold_;
  bool is_sort_;
  usint32_t K_;
  usint32_t M_;
  Evaluation* ground_truth_;
};

#endif /* EVALUATION_H_ */
