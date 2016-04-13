/*
 *    This file contains functions for evaluating different protein search tools.
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

#include "evaluation.hpp"
#include <algorithm>

void Evaluation::ShowReadResults(ofstream& fout) {
  for (ResultsType::const_iterator it = search_results.begin();
      it != search_results.end(); it++) {
    for (uint32_t i = 0; i < it->second.size(); i++) {
      fout << it->first << "\t" << it->second[i].protein_name << "\t"
          << it->second[i].identity << "\t" << it->second[i].aligned_len << "\t"
          << it->second[i].mismatch << "\t" << it->second[i].gap_open << "\t"
          << it->second[i].qs << "\t" << it->second[i].qe << "\t"
          << it->second[i].ps << "\t" << it->second[i].pe << "\t"
          << it->second[i].evalue << "\t" << it->second[i].bit_score << endl;
    }
  }
}

void Evaluation::ReadResults() {
  string strline;
  ifstream fin(results_file.c_str());
  cout << "read file: " << results_file << endl;
  while (!fin.eof()) {
    getline(fin, strline);
    if (strline.size() == 0)
      continue;
    if (strline[0] == '#')
      continue;
    vector < string > seg;
    string tmp;
    int tag = 0;
    for (uint32_t i = 0; i < strline.size(); i++) {
      if (strline[i] == '\t') {
        seg.push_back(tmp);
        tmp.clear();
        tag = 0;
      } else {
        if (strline[i] == ' ') {
          tag = 1;
        }
        if (tag == 0) {
          tmp += strline[i];
        }
      }
    }
    seg.push_back(tmp);

    double identity, evalue, bit_score;
    int aligned_len, mismatch, gap_open, qs, qe, ps, pe;
    sscanf(seg[2].c_str(), "%lf", &identity);
    aligned_len = atoi(seg[3].c_str());
    mismatch = atoi(seg[4].c_str());
    gap_open = atoi(seg[5].c_str());
    qs = atoi(seg[6].c_str());
    qe = atoi(seg[7].c_str());
    ps = atoi(seg[8].c_str());
    pe = atoi(seg[9].c_str());
    sscanf(seg[10].c_str(), "%lf", &evalue);
    sscanf(seg[11].c_str(), "%lf", &bit_score);

    if (evalue > evalue_threshold)
      continue;

    search_results[seg[0]].push_back(
        M8Results(seg[1], identity, aligned_len, mismatch, gap_open, qs, qe, ps,
                  pe, evalue, bit_score));
  }
  fin.close();
  if (is_sort) {
    cout << "sorting..." << endl;
    for (ResultsType::iterator it = search_results.begin();
        it != search_results.end(); it++) {
      sort(it->second.begin(), it->second.end(), M8Results::SORT_CMP_EValue);
    }
    ofstream fout(string(results_file + "_sort.txt").c_str());
    ShowReadResults(fout);
    fout.close();
    ////////////////////////////////
    //statistic for start position
    unordered_map < uint32_t, uint32_t > start_count;
    for (ResultsType::iterator it = search_results.begin();
        it != search_results.end(); ++it) {
      for (uint32_t i = 0; i < it->second.size() && i < 10; ++i) {
        start_count[it->second[0].qs]++;
      }
    }

    cout << string(results_file + "_start_count.txt") << endl;
    ofstream fsee(string(results_file + "_start_count.txt").c_str());
    for (unordered_map<uint32_t, uint32_t>::const_iterator it = start_count
        .begin(); it != start_count.end(); ++it) {
      fsee << it->first << " " << it->second << endl;
    }
    fsee.close();
  }
}

void Evaluation::CompareTopKHitProtein(ofstream& fout) {
  for (double evalue = 1e-64; evalue <= 100; evalue *= 100) {
    cout << "evalue = " << evalue << endl;
    uint32_t total = 0;
    uint32_t occur = 0;
    double hits = 0;
    for (ResultsType::const_iterator it = ground_truth->search_results.begin();
        it != ground_truth->search_results.end(); it++) {
      ResultsType::const_iterator it2 = search_results.find(it->first);
      for (uint32_t i = 0; i < it->second.size() && i < K; i++) {
        if (it->second[i].evalue > evalue)
          continue;
        total++;
        if (it2 == search_results.end())
          continue;

        for (uint32_t j = 0; j < it2->second.size() && j < M; j++) {
          if (it2->second[j].protein_name == it->second[i].protein_name) {
            occur++;
            break;
          }
        }
      }
    }
    cout << "total = " << total << endl;
    cout << "occur = " << occur << endl;
    cout << "accuracy = " << double(occur) / double(total) << endl;
    cout << "hits = " << hits << " " << double(hits) / (double) occur << endl;
    fout << double(occur) / double(total) << "\t";
  }
  fout << endl;
}

void Evaluation::CompareTopKHitProtein_INFO() {
  double evalue = 1e-30;
  cout << "evalue = " << evalue << endl;
  uint32_t total = 0;
  uint32_t occur = 0;
  double hits = 0;
  int tag = 0;

  ofstream fout("see_diff.txt");

  for (ResultsType::const_iterator it = ground_truth->search_results.begin();
      it != ground_truth->search_results.end(); it++) {
    ResultsType::const_iterator it2 = search_results.find(it->first);
    for (uint32_t i = 0; i < it->second.size() && i < K; i++) {
      if (it->second[i].evalue > evalue)
        continue;
      total++;
      if (it2 == search_results.end())
        continue;

      tag = 0;
      for (uint32_t j = 0; j < it2->second.size() && j < M; j++) {
        if (it2->second[j].protein_name == it->second[i].protein_name) {
          occur++;
          tag = 1;
          break;
        }
      }
      if (tag == 0) {
        fout << it->first << " " << it->second[i].protein_name << endl;
        DisplayResults(it->first.c_str(), "UNKNOW", it->second,
                       it->second.size(), 7, fout);
        fout << "---------------------------------------------------------"
            << endl;
        DisplayResults(it->first.c_str(), "UNKNOW", it2->second,
                       it2->second.size(), 7, fout);
        fout
            << "***************************************************************"
            << endl;
      }
    }
  }
  cout << "total = " << total << endl;
  cout << "occur = " << occur << endl;
  cout << "accuracy = " << double(occur) / double(total) << endl;
  cout << "hits = " << hits << " " << double(hits) / (double) occur << endl;
  fout << double(occur) / double(total) << "\t";
  fout.close();
}


