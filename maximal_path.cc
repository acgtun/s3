#include "maximal_path.h"
#include <limits.h>

namespace maximal_path {

void MaximalPath::TopologicalSortUtil(const usint32_t& v) {
  visited[v] = true;
  for (usint32_t i = 0; i < vLinkNodeSize[v]; i++)
    if (!visited[vLinkList[v][i]])
      TopologicalSortUtil(vLinkList[v][i]);

  top_order_.push_back(v);
}

void MaximalPath::TopologicalSort() {
  top_order_.clear();
  for (usint32_t i = 0; i < num_seg; i++)
    visited[i] = false;

  for (usint32_t i = 0; i < num_seg; i++)
    if (visited[i] == false)
      TopologicalSortUtil(i);

  //for (usint32_t i = 0; i < top_order_.size(); i++)
  // cout << i << " top_order " << top_order_[i] << endl;
}

int MaximalPath::WeightBetweenTwoHits(const CSegHit& a, const CSegHit& b) {
  usint32_t n = b.qs - a.qe - 1;
  usint32_t m = b.ps - a.pe - 1;
  if (m == n) {
    int score = 0;
    for (usint32_t i = 0; i < n; i++) {
      score += BLOSUM62[base[query_[a.qe + i] - 'A']][base[protein_[a.pe + i]
          - 'A']];
    }
    return score;
  }
  if (is_approximate_gaps) {
    return 0;
  } else {
    return global_alignment_.RunGlobalAlignment(&(query_[a.qe + 1]),
                                                &(protein_[a.pe + 1]),
                                                b.qs - a.qe - 1,
                                                b.ps - a.pe - 1);
  }
}

int MaximalPath::SegScore(const CSegHit& v) {
  int score = 0;
  for (usint32_t i = 0; i <= v.qe - v.qs; i++) {
    score += BLOSUM62[base[query_[v.qs + i] - 'A']][base[protein_[v.ps + i]
        - 'A']];
  }
  return score;
}

void MaximalPath::BuildGraph(const vector<CSegHit>& seg) {
  for (usint32_t i = 0; i < num_seg; i++) {
    vLinkNodeSize[i] = 0;
  }
  for (usint32_t i = 0; i < num_seg; i++) {
    for (usint32_t j = i + 1; j < num_seg; j++) {
      if (seg[j].ps < seg[i].pe || seg[j].qs < seg[i].qe) {
        continue;
      }
      if (seg[j].ps == seg[i].pe || seg[j].qs == seg[i].qe) {
        vLinkList[i][vLinkNodeSize[i]] = j;
        vEdgeWeight[i][vLinkNodeSize[i]] = gapopen_;
        vLinkNodeSize[i]++;
        continue;
      }
      if (max(seg[j].ps - seg[i].pe, seg[j].qs - seg[i].qe) < xdrop) {
        vLinkList[i][vLinkNodeSize[i]] = j;
        vEdgeWeight[i][vLinkNodeSize[i]] = WeightBetweenTwoHits(seg[i], seg[j]);
        vLinkNodeSize[i]++;
      }
    }
  }
  /*
   for (usint32_t i = 0; i < num_seg; i++) {
   cout << "i = " << i << ": ";
   for (usint32_t j = 0; j < vLinkNodeSize[i]; j++) {
   cout << "(" << vLinkList[i][j] << ", " << vEdgeWeight[i][j] << ")";
   }
   cout << endl;
   }*/
}

bool MaximalPath::FindSignificantPath(const char * query, const char * protein,
                                      const vector<CSegHit>& seg,
                                      M8Results& res) {
  //cout << "----------------START-------------------------" << endl;
  //cout << query << " " << endl;
  //cout << protein << endl;
  num_seg = seg.size();
  query_ = query;
  protein_ = protein;
  //for (usint32_t i = 0; i < num_seg; i++) {
  // cout << i << " " << seg[i].qs << " " << seg[i].qe << " " << seg[i].ps << " "
  // << seg[i].pe << endl;
  //}
  BuildGraph(seg);
  TopologicalSort();

  for (usint32_t i = 0; i < num_seg; i++) {
    seg_score_[i] = SegScore(seg[i]);
    score_[i] = seg_score_[i];
    //cout << i << " score " << seg_score_[i] << endl;
  }

  for (usint32_t i = 0; i < num_seg; i++) {
    parent_[i] = -1;
  }
  //cout << "Find path..." << endl;
  for (vector<usint32_t>::reverse_iterator rit = top_order_.rbegin();
      rit != top_order_.rend(); rit++) {
    usint32_t u = *rit, v;
    if (score_[u] != INT_MIN) {
      for (usint32_t i = 0; i < vLinkNodeSize[u]; i++) {
        v = vLinkList[u][i];
        if (score_[v] < score_[u] + vEdgeWeight[u][v] + seg_score_[v]) {
          score_[v] = score_[u] + vEdgeWeight[u][v] + seg_score_[v];
          parent_[v] = u;
        }
      }
    }
  }

  int max_score = -1;
  int node = 0;
  for (usint32_t i = 0; i < num_seg; i++) {
    //cout << i << "score_f " << score_[i] << endl;
    if (score_[i] > max_score) {
      node = i;
      max_score = score_[i];
    }
  }
  //cout << "max_node = " << node << endl;
  // for (usint32_t i = 0; i < num_seg; i++) {
  // cout << i << " parent " << parent_[i] << endl;
  //}
  int j = node;
  int len = seg[j].qe - seg[j].qs + 1;
  while (parent_[j] != -1) {
    len += max(seg[j].qs - seg[parent_[j]].qe + 1,
               seg[j].ps - seg[parent_[j]].pe + 1);
    j = parent_[j];
    len += (seg[j].qe - seg[j].qs + 1);
  }

  if (score_[node] >= evalue_->score_by_evalue_[1]) {
    //cout << "yes " << seg[j].qs << " " << seg[node].qe << " " << seg[j].ps
    //  << " " << seg[j].pe << endl;

    res = M8Results("", 1, seg[0].qe - seg[0].qs + 1, 0, 0, seg[j].qs,
                    seg[node].qe, seg[j].ps, seg[node].pe,
                    evalue_->GetEvalue(score_[node], 1),
                    evalue_->GetBitScore(score_[node], 1));
    return true;
  }
  return false;
  //cout << "------------------------END-----------------" << endl;
}

}  //namespace maximal_path

