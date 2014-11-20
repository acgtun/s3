#ifndef QUERY_SEARCH_H_
#define QUERY_SEARCH_H_

#include "./util/bio_util.h"
#include "./util/option.h"
#include "./index/index.h"
#include "hits_matrix.h"
#include "./util/evalue.h"
#include "./seqalign/local_alignment.h"

using local_alignment::LocalAlignment;
using hits_matrix::HitsMatrix;

class QuerySearch {
 public:
  QuerySearch(const Index* index, const Evalue* evalue,
              const usint32_t& max_query_length)
      : hits_matrix_(evalue, max_query_length,
                     index->database_.max_sequence_length_),
        index_(index),
        evalue_(evalue) {
    protein_hits_threshold_ = 0;
    query_name_ = NULL;

    Option::GetOption("-tpro", num_of_candidate_protein, 100);

    fres.open("hfresults.out");
    fseeprotein.open("protein_see.txt");
  }

  ~QuerySearch() {
    fres.close();
    fseeprotein.close();
  }
  void Search(const char* query_seq, const char* query_name);

 private:
  void GetProteinScoreThreshold();
  void ShowHits();
  void ShowProteinHasHits();
  ofstream fres;
  ofstream fseeprotein;
  const char* query_name_;

  map<usint32_t, vector<CSegHit> > protein_hits;
  usint32_t protein_hits_threshold_;
  usint32_t num_of_candidate_protein;
  priority_queue<usint32_t, vector<usint32_t>, greater<usint32_t> > top_protein_;

  vector<M8Results> aligned_results;

  HitsMatrix hits_matrix_;
  const Index *index_;
  const Evalue* evalue_;
};

#endif /* QUERY_SEARCH_H_ */
