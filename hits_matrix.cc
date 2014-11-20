#include "hits_matrix.h"
#include <limits.h>

namespace hits_matrix {

void HitsMatrix::ShowHits(map<int, vector<CSegHit> >& seg) {
  cout << "------------------------------------------" << endl;
  for (map<int, vector<CSegHit> >::iterator it = seg.begin(); it != seg.end();
      it++) {
    cout << "diag " << it->first << endl;
    for (usint32_t i = 0; i < it->second.size(); i++) {
      cout << "(" << it->second[i].qs << " " << it->second[i].qe << " "
          << it->second[i].ps << " " << it->second[i].pe << ")";
    }
    cout << endl;
  }
}

bool SortByQueryStartPosition(const CSegHit& a, const CSegHit& b) {
  return a.qs < b.qs;
}

int HitsMatrix::DiagID(const usint32_t & qs, const usint32_t & ps) {
  return (int) qs - (int) ps;
}

void HitsMatrix::SplitToDiag(const vector<CSegHit>& hits,
                             map<int, vector<CSegHit> >& seg) {
  seg.clear();
  for (usint32_t i = 0; i < hits.size(); i++) {
    seg[DiagID(hits[i].qs, hits[i].ps)].push_back(hits[i]);
  }
}

void HitsMatrix::ChainFilter(map<int, vector<CSegHit> >& seg,
                             map<int, vector<CSegHit> >& seg_new) {
  seg_new.clear();
  for (map<int, vector<CSegHit> >::iterator it = seg.begin(); it != seg.end();
      it++) {
    sort(it->second.begin(), it->second.end(), SortByQueryStartPosition);
    CSegHit one_seg(it->second[0]);
    for (usint32_t j = 1; j < it->second.size(); j++) {
      if (it->second[j].qs <= one_seg.qe + 5) {
        one_seg.qe = it->second[j].qe;
        one_seg.pe = it->second[j].pe;
      } else {
        seg_new[it->first].push_back(one_seg);
        one_seg = CSegHit(it->second[j]);
      }
    }
    seg_new[it->first].push_back(one_seg);
  }
}

void HitsMatrix::UngappedExtension(map<int, vector<CSegHit> >& seg) {
  for (map<int, vector<CSegHit> >::iterator it = seg.begin(); it != seg.end();
      it++) {
    for (usint32_t i = 0; i < it->second.size(); i++) {
      usint32_t qe = it->second[i].qe + 1;
      usint32_t pe = it->second[i].pe + 1;
      int scoreMax = SegScore(it->second[i]);
      int score = scoreMax;
      while (qe < strlen(query_) && pe < strlen(protein_)) {
        score += BLOSUM62[base[query_[qe] - 'A']][base[protein_[pe] - 'A']];
        if (score < scoreMax - D) {
          it->second[i].qe = qe - 1;
          it->second[i].pe = pe - 1;
          break;
        }
        qe++;
        pe++;
      }

      usint32_t qs = it->second[i].qs, ps = it->second[i].ps;
      score = scoreMax;
      while (qs > 0 && ps > 0) {
        qs--;
        ps--;
        score += BLOSUM62[base[query_[qs] - 'A']][base[protein_[ps] - 'A']];
        if (score < scoreMax - D) {
          it->second[i].qs = qs + 1;
          it->second[i].ps = ps + 1;
          break;
        }
      }
    }
  }
}

int HitsMatrix::SegScore(const CSegHit& v) {
  int score = 0;
  for (usint32_t i = 0; i <= v.qe - v.qs; i++) {
    score += BLOSUM62[base[query_[v.qs + i] - 'A']][base[protein_[v.ps + i]
        - 'A']];
  }
  return score;
}

void HitsMatrix::SelectDiags(map<int, vector<CSegHit> >& seg,
                             vector<M8Results>& aligned_results) {
  int diagEnd = seg.begin()->first;
  vector<CSegHit> seg_ptr;
  M8Results res;
  for (map<int, vector<CSegHit> >::iterator it = seg.begin(); it != seg.end();
      it++) {
    if (it->first - diagEnd < diags_diff) {
      for (usint32_t i = 0; i < it->second.size(); i++) {
        seg_ptr.push_back(it->second[i]);
      }
      diagEnd = it->first;
    } else {
      if (using_local_alignment) {
        usint32_t x1 = INT_MAX;
        usint32_t x2 = 0;
        usint32_t y1 = INT_MAX;
        usint32_t y2 = 0;
        for (usint32_t i = 0; i < seg_ptr.size(); i++) {
          if (seg_ptr[i].qs < x1)
            x1 = seg_ptr[i].qs;
          if (seg_ptr[i].qe > x2)
            x2 = seg_ptr[i].qe;

          if (seg_ptr[i].ps < y1)
            y1 = seg_ptr[i].ps;
          if (seg_ptr[i].pe > y2)
            y2 = seg_ptr[i].pe;
        }
        if (local_alignment_.RunLocalAlignment(&(query_[x1]), &(protein_[y1]),
                                               x2 - x1 + 1, y2 - y1 + 1, res)) {
          res.protein_name = protein_name_;
          aligned_results.push_back(res);
        }
      } else {
        sort(seg_ptr.begin(), seg_ptr.end(), SortByQueryStartPosition);
        if (maximal_path_.FindSignificantPath(query_, protein_, seg_ptr, res)) {
          res.protein_name = protein_name_;
          aligned_results.push_back(res);
        }
      }
      //FindPath(seg_ptr);
      ////////////////
      seg_ptr.clear();
      diagEnd = it->first;
      for (usint32_t i = 0; i < it->second.size(); i++) {
        seg_ptr.push_back(it->second[i]);
      }
    }
  }
  if (seg_ptr.size() != 0) {

    //cout << "opt_->using_local_alignment = " << opt_->using_local_alignment
    //   << endl;
    if (using_local_alignment) {
      usint32_t x1 = INT_MAX;
      usint32_t x2 = 0;
      usint32_t y1 = INT_MAX;
      usint32_t y2 = 0;
      for (usint32_t i = 0; i < seg_ptr.size(); i++) {
        if (seg_ptr[i].qs < x1)
          x1 = seg_ptr[i].qs;
        if (seg_ptr[i].qe > x2)
          x2 = seg_ptr[i].qe;

        if (seg_ptr[i].ps < y1)
          y1 = seg_ptr[i].ps;
        if (seg_ptr[i].pe > y2)
          y2 = seg_ptr[i].pe;
      }
      //  cout << x1 << endl;
      // cout << x2 << endl;
      //cout << y1 << endl;
      // cout << y2 << endl;
      if (local_alignment_.RunLocalAlignment(&(query_[x1]), &(protein_[y1]),
                                             x2 - x1 + 1, y2 - y1 + 1, res)) {
        res.protein_name = protein_name_;
        aligned_results.push_back(res);
      }
    } else {
      sort(seg_ptr.begin(), seg_ptr.end(), SortByQueryStartPosition);
      if (maximal_path_.FindSignificantPath(query_, protein_, seg_ptr, res)) {
        res.protein_name = protein_name_;
        aligned_results.push_back(res);
      }
    }
  }
}

void HitsMatrix::AnalyzeMatrix(const char* query, const char* protein,
                               const string& protein_name,
                               const vector<CSegHit>& hits,
                               vector<M8Results>& aligned_results) {
  query_ = query;
  protein_ = protein;
  protein_name_ = protein_name;
  map<int, vector<CSegHit> > seg;
  map<int, vector<CSegHit> > seg_new;

  SplitToDiag(hits, seg);
  ChainFilter(seg, seg_new);
  UngappedExtension(seg_new);
  ChainFilter(seg_new, seg);
  SelectDiags(seg, aligned_results);
}

}  //namespace hits_matrix

