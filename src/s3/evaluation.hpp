/*
 *    This is the header file for evaluation.
 *
 *    Copyright (C) 2015 University of Southern California
 *
 *    Authors: Haifeng Chen and Ting Chen
 *
 *    This file is part of S3.
 *
 *    S3 is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    S3 is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with S3.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef EVALUATION_H_
#define EVALUATION_H_

#include "util.hpp"

typedef unordered_map<string, vector<M8Results> > ResultsType;

class Evaluation {
 public:
  Evaluation(const string& _results_file, const Evaluation* _ground_truth,
             const double& _evalue_threshold, const bool& _is_sort,
             const int& _K, const int& _M)
      : results_file(_results_file),
        ground_truth(_ground_truth),
        evalue_threshold(_evalue_threshold),
        is_sort(_is_sort),
        K(_K),
        M(_M) {
  }
  ~Evaluation() {
  }

  void ReadResults();
  void ShowReadResults(ofstream& fout);
  void CompareTopKHitProtein(ofstream& fout);
  void CompareTopKHitProtein_INFO();

  string results_file;
  ResultsType search_results;
  const Evaluation* ground_truth;

  double evalue_threshold;
  bool is_sort;
  uint32_t K;
  uint32_t M;
};

#endif /* EVALUATION_H_ */
