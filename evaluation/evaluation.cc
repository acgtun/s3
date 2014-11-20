#include "evaluation.h"
#include <algorithm>

void Evaluation::ShowReadResults(ofstream& fout) {
  for (map<string, vector<M8Results> >::iterator it = search_results_.begin();
       it != search_results_.end(); it++) {
    for (usint32_t i = 0; i < it->second.size(); i++) {
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
  ifstream fin(results_file_.c_str());
  cout << "read file: " << results_file_ << endl;
  while (!fin.eof()) {
    getline(fin, strline);
    if (strline.size() == 0)
      continue;
    if (strline[0] == '#')
      continue;
    vector < string > seg;
    string tmp;
    int tag = 0;
    for (usint32_t i = 0; i < strline.size(); i++) {
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

    if (evalue > evalue_threshold_)
      continue;

    search_results_[seg[0]].push_back(
        M8Results(seg[1], identity, aligned_len, mismatch, gap_open, qs, qe, ps,
                  pe, evalue, bit_score));
  }
  fin.close();
  if (is_sort_) {
    cout << "sorting..." << endl;
    for (map<string, vector<M8Results> >::iterator it = search_results_.begin();
         it != search_results_.end(); it++) {
      sort(it->second.begin(), it->second.end(), M8Results::SORT_CMP_EValue);
    }
    ofstream fout(string(results_file_ + "_sort.txt").c_str());
    ShowReadResults(fout);
    fout.close();
  }
}

void Evaluation::CompareTopKHitProtein(ofstream& fout) {
  for (double evalue = 1e-64; evalue <= 100; evalue *= 100) {
    cout << "evalue = " << evalue << endl;
    usint32_t total = 0;
    usint32_t occur = 0;
    double hits = 0;
    for (map<string, vector<M8Results> >::iterator it = ground_truth_
         ->search_results_.begin(); it != ground_truth_->search_results_.end();
         it++) {
      map<string, vector<M8Results> >::iterator it2 = search_results_.find(
          it->first);
      for (usint32_t i = 0; i < it->second.size() && i < K_; i++) {
        if (it->second[i].evalue > evalue)
          continue;
        total++;
        if (it2 == search_results_.end())
          continue;

        for (usint32_t j = 0; j < it2->second.size() && j < M_; j++) {
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
  usint32_t total = 0;
  usint32_t occur = 0;
  double hits = 0;
  int tag = 0;

  ofstream fout("see_diff.txt");

  for (map<string, vector<M8Results> >::iterator it = ground_truth_
       ->search_results_.begin(); it != ground_truth_->search_results_.end();
       it++) {
    map<string, vector<M8Results> >::iterator it2 = search_results_.find(
        it->first);
    for (usint32_t i = 0; i < it->second.size() && i < K_; i++) {
      if (it->second[i].evalue > evalue)
        continue;
      total++;
      if (it2 == search_results_.end())
        continue;

      tag = 0;
      for (usint32_t j = 0; j < it2->second.size() && j < M_; j++) {
        if (it2->second[j].protein_name == it->second[i].protein_name) {
          occur++;
          tag = 1;
          break;
        }
      }
      if(tag == 0) {
        fout << it->first << " " << it->second[i].protein_name << endl;      
        DisplayResults(it->first.c_str(), "UNKNOW", it->second, 7, fout);
        fout << "---------------------------------------------------------" << endl;
        DisplayResults(it->first.c_str(), "UNKNOW", it2->second, 7, fout);
        fout << "***************************************************************" << endl;
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
