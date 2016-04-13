/*
 *    This is the header file local alignment.
 *
 *    Copyright (C) 2016 University of Southern California
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

#ifndef LOCALALIGNMENT_H_
#define LOCALALIGNMENT_H_

#include "util.hpp"
#include "evalue.hpp"

namespace local_alignment {

using std::max;

#define DIAG ('a')
#define UP ('b')
#define LEFT ('c')
#define STOPorSTART ('d')

class LocalAlignment {
 public:
  LocalAlignment(Evalue* _evalue, const uint32_t& _max_rows,
                 const uint32_t& _max_cols, const int& _gapopen,
                 const int& _gapextension);
  ~LocalAlignment();

  bool RunLocalAlignment(const string& U, const string& V,
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

  Evalue* evalue;

  uint32_t max_rows;
  uint32_t max_cols;

  int gapopen;
  int gapextension;

  void DisplayAlignment();
  void stringReverse(vector<char>& str, const int & n);
  int MaxOfFour(const int & s1, const int & s2, const int & s3, const int & s4);
  char Direction(const int & s1, const int & s2, const int & s3,
                 const int & s4);
};

}  // namespace local_alignment

#endif /* LOCALALIGNMENT_H_ */
