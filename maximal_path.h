#ifndef MAXIMAL_PATH_H_
#define MAXIMAL_PATH_H_

#include "./seqalign/global_alignment.h"
#include "./util/option.h"
#include "./index/index.h"
#include "./util/fasta_file.h"
#include <stack>
#include <list>
#include "./util/bio_util.h"
#include "./util/evalue.h"

using global_alignment::GlobalAlignment;

namespace maximal_path {

#define MAX_ALIGNMENT_LEN 6000

struct CSegHit {
  CSegHit(const usint32_t & qstart, const usint32_t & qend,
          const usint32_t & pstart, const usint32_t& pend)
      : qs(qstart),
        qe(qend),
        ps(pstart),
        pe(pend) {
  }
  CSegHit(const CSegHit& seg)
      : qs(seg.qs),
        qe(seg.qe),
        ps(seg.ps),
        pe(seg.pe) {
  }

  usint32_t qs;
  usint32_t qe;
  usint32_t ps;
  usint32_t pe;
  //int score;
};

class MaximalPath {
 public:
  MaximalPath(const Evalue* evalue)
      : global_alignment_(),
        evalue_(evalue) {
    Option::GetOption("-fast", is_approximate_gaps, 0);
    Option::GetOption("-x", xdrop, 50);
    Option::GetOption("-gopen", gapopen_, 50);

    query_ = NULL;
    protein_ = NULL;

    vLinkList.resize(MAX_ALIGNMENT_LEN);
    vLinkNodeSize.resize(MAX_ALIGNMENT_LEN);
    vEdgeWeight.resize(MAX_ALIGNMENT_LEN);
    seg_score_.resize(MAX_ALIGNMENT_LEN);

    for (usint32_t i = 0; i < MAX_ALIGNMENT_LEN; i++) {
      vLinkList[i].resize(MAX_ALIGNMENT_LEN);
      vEdgeWeight[i].resize(MAX_ALIGNMENT_LEN);
    }

    visited.resize(MAX_ALIGNMENT_LEN);
    parent_.resize(MAX_ALIGNMENT_LEN);
    score_.resize(MAX_ALIGNMENT_LEN);
  }

  ~MaximalPath() {
  }
  bool FindSignificantPath(const char * query, const char * protein,
                           const vector<CSegHit>& seg, M8Results& res);

 private:
  void BuildGraph(const vector<CSegHit>& seg);
  void TopologicalSort();
  void TopologicalSortUtil(const usint32_t& v);
  int WeightBetweenTwoHits(const CSegHit& a, const CSegHit& b);
  int SegScore(const CSegHit& v);

  vector<vector<usint32_t> > vLinkList;
  vector<usint32_t> vLinkNodeSize;
  vector<vector<int> > vEdgeWeight;
  vector<int> seg_score_;

  vector<usint32_t> top_order_;
  vector<bool> visited;

  vector<int> parent_;

  vector<int> score_;

  const char* query_;
  const char* protein_;

  usint32_t num_seg;

  GlobalAlignment global_alignment_;
  const Evalue* evalue_;

  int is_approximate_gaps;
  usint32_t xdrop;
  int gapopen_;

};

}
#endif /* MAXIMAL_PATH_H_ */
