#include "query_search.h"
#include <stdlib.h>
#include <vector>

void QuerySearch::GetProteinScoreThreshold() {
  for (map<usint32_t, vector<CSegHit> >::iterator it = protein_hits.begin();
      it != protein_hits.end(); it++) {
    if (it->second.size() >= 1) {
      top_protein_.push(it->second.size());
      if (top_protein_.size() == num_of_candidate_protein + 1) {
        top_protein_.pop();
      }
    }
  }
  if (top_protein_.size() <= num_of_candidate_protein) {
    protein_hits_threshold_ = 0;
  } else {
    protein_hits_threshold_ = top_protein_.top();
  }

  while (!top_protein_.empty()) {
    top_protein_.pop();
  }
}

void QuerySearch::Search(const char* query_seq, const char* query_name) {
  query_name_ = query_name;
  protein_hits.clear();
  aligned_results.clear();
  usint32_t numSeed = strlen(query_seq) - HASHAALEN + 1;
  usint32_t hashValue, Kmer, S, E;

  for (usint32_t j = 0; j < numSeed; j++) {
    /* (1) Get Hash Value of the Kmer (here Kmer is seed) */
    hashValue = GetHashValue(&(query_seq[j]));

    /* (2) Get Nearest Kmers from the kmer_nearest_ table */
    for (usint64_t c = index_->kmer_nearest_.counter_[hashValue];
        c < index_->kmer_nearest_.counter_[hashValue + 1]; c++) {
      /* (3) Get the Kmer */
      Kmer = index_->kmer_nearest_.index_[c];
      //score = index_->kmer_nearest_.score_[c];

      /* (3) Get all the protein database position which occurs the Kmer */
      S = index_->kmer_pos_.counter_[Kmer];
      E = index_->kmer_pos_.counter_[Kmer + 1];

      for (usint64_t c2 = S; c2 < E; c2++) {
        protein_hits[index_->kmer_pos_.index_pro_id_[c2]].push_back(
            CSegHit(j, j + HASHAALEN - 1, index_->kmer_pos_.index_pro_pos_[c2],
                    index_->kmer_pos_.index_pro_pos_[c2] + HASHAALEN - 1));
      }
    }
  }

  /* Get Protein Score Threshold */
  GetProteinScoreThreshold();

  for (map<usint32_t, vector<CSegHit> >::iterator it = protein_hits.begin();
      it != protein_hits.end(); it++) {
    if (it->second.size() < protein_hits_threshold_)
      continue;
    hits_matrix_.AnalyzeMatrix(query_seq,
                               index_->database_.sequences_[it->first],
                               index_->database_.sequences_names_[it->first],
                               it->second, aligned_results);
  }

  DisplayResults(query_name, index_->database_.GetFilePath().c_str(),
                 aligned_results, 7, fres);
}
