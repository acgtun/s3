#include "k_mer_nearest.h"

KmerNearest::KmerNearest(const KmerPosition* kmer_pos)
    : kmer_pos_(kmer_pos) {
  counter_ = NULL;
  index_ = NULL;
  size_of_counter_ = 0;
  size_of_index_ = 0;
  kmer_score_threshold_ = 30;
}

KmerNearest::~KmerNearest() {
}

void KmerNearest::BuildNearestKmer() {
  INFO("Build Nearest Kmer...");
  MEMORY_ALLOCATE_CHECK(
      index_ = (usint32_t * ) malloc(sizeof(usint32_t) * size_of_index_));
  //MEMORY_ALLOCATE_CHECK(score_ = (int * ) malloc(sizeof(int) * size_of_index_));
  usint32_t hashValueCand = 0;
  char querySeq[HASHAALEN + 1], strCandPep[HASHAALEN + 1];
  for (usint32_t i = 0; i < kmer_pos_->size_of_counter_; i++) {
    if (kmer_pos_->counter_[i] == kmer_pos_->counter_[i + 1])
      continue;
    DeCodeAA(i, querySeq);
    k_longest_path_.UpdateWeight(querySeq);
    bool bPath = true;
    int pathID = 0, rawScore = 0;
#ifdef TESTCODE
    set < string > find_candidate_test;
    set<usint32_t> hashValueCand_test;
#endif
    while (bPath && pathID < MAX_K) {
      bPath = k_longest_path_.FindCandidate(strCandPep, pathID, rawScore);
#ifdef TESTCODE
      if (find_candidate_test.find(strCandPep) != find_candidate_test.end()) {
        ERROR_INFO("The Find k Longest Path has error");
      } else {
        find_candidate_test.insert(strCandPep);
      }
#endif
      if (pathID == 0 && rawScore < kmer_score_threshold_) {
        if (strcmp(strCandPep, querySeq) != 0) {
          ERROR_INFO("The Find k Longest Path has error");
        }
        index_[counter_[i]] = i;
        //score_[counter_[i]] = rawScore;
        counter_[i]++;
        break;
      }
      if (rawScore < kmer_score_threshold_)
        break;
      hashValueCand = GetHashValue(strCandPep);
#ifdef TESTCODE
      if (hashValueCand_test.find(hashValueCand) != hashValueCand_test.end()) {
        ERROR_INFO("The hash function has error~");
      } else
        hashValueCand_test.insert(hashValueCand);
#endif
      index_[counter_[hashValueCand]] = i;
     // score_[counter_[hashValueCand]] = rawScore;
      counter_[hashValueCand]++;
      pathID++;
    }
  }

  for (usint32_t i = size_of_counter_ - 1; i >= 1; i--) {
    counter_[i] = counter_[i - 1];
  }
  counter_[0] = 0;
}

void KmerNearest::CountNearestKmer() {
  INFO("Count Nearest Kmer...");
  size_of_counter_ = GetPower(20, HASHAALEN);
  INFO("The size of Kmer Nearest Counter is ", size_of_counter_);
  MEMORY_ALLOCATE_CHECK(
      counter_ = (usint64_t * ) malloc(
          sizeof(usint64_t) * (size_of_counter_ + 1)));
  memset(counter_, 0, sizeof(usint64_t) * (size_of_counter_ + 1));

  usint32_t hashValueCand = 0;
  char querySeq[HASHAALEN + 1], strCandPep[HASHAALEN + 1];
  for (usint32_t i = 0; i < kmer_pos_->size_of_counter_; i++) {
    if (kmer_pos_->counter_[i] == kmer_pos_->counter_[i + 1])
      continue;
    DeCodeAA(i, querySeq);
    k_longest_path_.UpdateWeight(querySeq);
    bool bPath = true;
    int pathID = 0, rawScore = 0;
#ifdef TESTCODE
    set < string > find_candidate_test;
    set<usint32_t> hashValueCand_test;
#endif
    while (bPath && pathID < MAX_K) {
      bPath = k_longest_path_.FindCandidate(strCandPep, pathID, rawScore);
#ifdef TESTCODE
      if (find_candidate_test.find(strCandPep) != find_candidate_test.end()) {
        ERROR_INFO("The Find k Longest Path has error");
      } else {
        find_candidate_test.insert(strCandPep);
      }
#endif
      if (pathID == 0 && rawScore < kmer_score_threshold_) {
        if (strcmp(strCandPep, querySeq) != 0) {
          printf("CAND=%s\n", strCandPep);
          printf("QUERY=%s\n", querySeq);
          cout << "i = " << i << " " << "rawScore = " << rawScore << endl;
          ERROR_INFO("The Find k Longest Path has error");
        }
        counter_[i]++;
      //  cout << i << endl;
        break;
      }
      if (rawScore < kmer_score_threshold_)
        break;
      hashValueCand = GetHashValue(strCandPep);
#ifdef TESTCODE
      if (hashValueCand_test.find(hashValueCand) != hashValueCand_test.end()) {
        ERROR_INFO("The hash function has error~");
      } else
        hashValueCand_test.insert(hashValueCand);
#endif
      counter_[hashValueCand]++;
      pathID++;
    }
  }

  // 2^32 / 20^6
  // Here we should be careful about the sum since the sum may more
  // than 2^32 if the average nubmer of nearest for each query is
  // more than certain number.
  //cout << "size_of_counter_ = " << size_of_counter_ << endl;
  for (usint32_t i = 1; i <= size_of_counter_; i++) {
    //cout << i << " " << counter_[i]  << " ";
    counter_[i] += counter_[i - 1];
    //cout << counter_[i]  << endl;
  }
  //cout << "finish.." << endl;
  size_of_index_ = counter_[size_of_counter_];
  INFO("The size of Kmer Nearest Index is ", size_of_index_);

  for (usint32_t i = size_of_counter_ - 1; i >= 1; i--) {
    counter_[i] = counter_[i - 1];
  }
  counter_[0] = 0;
}
