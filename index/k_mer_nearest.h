#ifndef KMERNEAREST_H_
#define KMERNEAREST_H_

#include "./../util/bio_util.h"
#include "k_mer_position.h"
#include "k_longest_path.h"

using longest_path::KLongestPath;

class KmerNearest {
 public:
  KmerNearest(const KmerPosition* kmer_pos);
  ~KmerNearest();

  void set_kmer_score_threshold_(const int& kmer_score_threshold) {
    kmer_score_threshold_ = kmer_score_threshold;
  }

  void BuildNearestKmer();
  void CountNearestKmer();

 public:
  usint64_t* counter_;
  usint32_t* index_;
 // int* score_;
  usint32_t size_of_counter_;
  usint64_t size_of_index_;

 private:
  const KmerPosition* kmer_pos_;
  KLongestPath k_longest_path_;
  int kmer_score_threshold_;
};

#endif /* KMERNEAREST_H_ */
