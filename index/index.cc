#include "index.h"
#include <stdlib.h>
#include <string>
#include <stdio.h>
#include <string.h>

void Index::GetThresholdScore() {
  INFO("Get Threshold Score...");
  map<int, usint64_t> org_count;
  for (usint32_t i = 0; i < AA20.size(); i++) {
    for (usint32_t j = 0; j < AA20.size(); j++) {
      org_count[BLOSUM62[i][j]]++;
    }
  }

  map<int, usint64_t> new_count;
  map<int, usint64_t> old_count(org_count);
  for (int i = 1; i < HASHAALEN; i++) {
    new_count.clear();
    for (map<int, usint64_t>::iterator it1 = old_count.begin();
        it1 != old_count.end(); it1++) {
      for (map<int, usint64_t>::iterator it2 = org_count.begin();
          it2 != org_count.end(); it2++) {
        new_count[it1->first + it2->first] += it1->second * it2->second;
      }
    }
    old_count.clear();
    old_count = new_count;
  }

  usint64_t sum = 0;
  INFO("Calculate Threshold Score...");
  for (map<int, usint64_t>::iterator it = old_count.begin();
      it != old_count.end(); it++) {
    //cout << it->first << " " << it->second << endl;
    sum += it->second;
  }

  //cout << database_.num_of_characters_ << endl;
  double p = (double) expected_hits_for_query_
      / (double) database_.num_of_characters_;
  INFO("The probability for each pair that has score larger than T is ", p);
  usint64_t threshold = 0;
  double kmer_score_threshold = 0;
  for (map<int, usint64_t>::reverse_iterator it = old_count.rbegin();
      it != old_count.rend(); it++) {
    //cout << it->first << " " << it->second << endl;
    threshold += it->second;
    if ((double) threshold / (double) sum > p) {
      kmer_score_threshold = it->first;
      break;
    }
  }
  if (kmer_score_threshold < 0)
    kmer_score_threshold = 0;
  INFO("The threshold score is", kmer_score_threshold);
  kmer_nearest_.set_kmer_score_threshold_(kmer_score_threshold);
}

void Index::BuildIndex() {
  GetThresholdScore();

  TIME_INFO(kmer_pos_.CountKmerPosition(), "Count Kmer Position");
  TIME_INFO(kmer_pos_.BuildKmerPosition(), "Build Kmer Position");

  TIME_INFO(kmer_nearest_.CountNearestKmer(), "Count Neareast Kmer");
  TIME_INFO(kmer_nearest_.BuildNearestKmer(), "Build Nearest Kmer");

  TIME_INFO(TestIndex(), "Write Test Index");
  TIME_INFO(WriteIndex(), "Write Index");
}

void Index::WriteIndex() {
  FILE * fout = fopen(index_file_.c_str(), "wb");
  INFO("Write index to", index_file_.c_str());

  fwrite(&(database_.num_of_characters_), sizeof(usint64_t), 1, fout);
  fwrite(&(database_.num_of_sequences_), sizeof(usint32_t), 1, fout);
  fwrite(&(database_.max_sequence_length_), sizeof(usint32_t), 1, fout);
  for (usint32_t i = 0; i < database_.num_of_sequences_; i++) {
    usint32_t l = database_.sequences_names_[i].size();
    fwrite(&l, sizeof(usint32_t), 1, fout);
    fwrite(&(database_.sequences_names_[i][0]), sizeof(char), l, fout);
    l = strlen(database_.sequences_[i]);
    fwrite(&l, sizeof(usint32_t), 1, fout);
    fwrite(database_.sequences_[i], sizeof(char), l, fout);
  }

  fwrite(&(kmer_nearest_.size_of_counter_), sizeof(usint32_t), 1, fout);
  fwrite(&(kmer_nearest_.size_of_index_), sizeof(usint64_t), 1, fout);
  fwrite(kmer_nearest_.counter_, sizeof(usint64_t),
         kmer_nearest_.size_of_counter_ + 1, fout);
  fwrite(kmer_nearest_.index_, sizeof(usint32_t), kmer_nearest_.size_of_index_,
         fout);
  // fwrite(kmer_nearest_.score_, sizeof(int), kmer_nearest_.size_of_index_, fout);

  fwrite(&(kmer_pos_.size_of_counter_), sizeof(usint32_t), 1, fout);
  fwrite(&(kmer_pos_.size_of_index_), sizeof(usint64_t), 1, fout);
  fwrite(kmer_pos_.counter_, sizeof(usint64_t), kmer_pos_.size_of_counter_ + 1,
         fout);
  fwrite(kmer_pos_.index_pro_id_, sizeof(usint32_t), kmer_pos_.size_of_index_,
         fout);
  fwrite(kmer_pos_.index_pro_pos_, sizeof(usint32_t), kmer_pos_.size_of_index_,
         fout);

  fclose(fout);
}

void Index::TestIndex() {
  INFO("Write Test File...");
  ofstream fp("table_position.txt");
  char querySeq[HASHAALEN + 1];
  for (usint32_t i = 0; i < kmer_pos_.size_of_counter_; i++) {
    if (kmer_pos_.counter_[i] == kmer_pos_.counter_[i + 1])
      continue;
    DeCodeAA(i, querySeq);
    fp << i << " " << querySeq << " " << kmer_pos_.counter_[i] << endl;
    for (usint64_t j = kmer_pos_.counter_[i]; j < kmer_pos_.counter_[i + 1];
        j++) {
      fp << database_.sequences_names_[kmer_pos_.index_pro_id_[j]] << " "
          << kmer_pos_.index_pro_pos_[j] << " ";
      fp
          << database_.sequences_[kmer_pos_.index_pro_id_[j]][kmer_pos_
              .index_pro_pos_[j]]
          << database_.sequences_[kmer_pos_.index_pro_id_[j]][kmer_pos_
              .index_pro_pos_[j] + 1]
          << database_.sequences_[kmer_pos_.index_pro_id_[j]][kmer_pos_
              .index_pro_pos_[j] + 2]
          << database_.sequences_[kmer_pos_.index_pro_id_[j]][kmer_pos_
              .index_pro_pos_[j] + 3]
          << database_.sequences_[kmer_pos_.index_pro_id_[j]][kmer_pos_
              .index_pro_pos_[j] + 4]
          << database_.sequences_[kmer_pos_.index_pro_id_[j]][kmer_pos_
              .index_pro_pos_[j] + 5] << endl;
    }
  }
  fp.close();

  ofstream fn("table_nearest.txt");
  char queryCand[HASHAALEN + 1];
  for (usint32_t i = 0; i < kmer_nearest_.size_of_counter_; i++) {
    if (kmer_nearest_.counter_[i] == kmer_nearest_.counter_[i + 1])
      continue;
    DeCodeAA(i, querySeq);
    fn << querySeq << " " << i << " (";
    for (usint64_t j = kmer_nearest_.counter_[i];
        j < kmer_nearest_.counter_[i + 1]; j++) {
      DeCodeAA(kmer_nearest_.index_[j], queryCand);
      fn << queryCand << " " << kmer_nearest_.index_[j] << " " << "SCORE"/*kmer_nearest_.score_[j]*/
      << " ,";
    }
    fn << ")" << endl;
  }
  fn.close();
}

void Index::ReadIndex() {
  INFO("Read index from", index_file_.c_str());
  FILE * fin = fopen(index_file_.c_str(), "rb");
  FILE_OPEN_CHECK(fin);

  FREAD_CHECK(fread(&(database_.num_of_characters_), sizeof(usint64_t), 1, fin),
              1);
  FREAD_CHECK(fread(&(database_.num_of_sequences_), sizeof(usint32_t), 1, fin),
              1);
  FREAD_CHECK(
      fread(&(database_.max_sequence_length_), sizeof(usint32_t), 1, fin), 1);

  database_.sequences_names_.resize(database_.num_of_sequences_);

  char proName[MAX_LINE_LEN];
  usint32_t l;
  MEMORY_ALLOCATE_CHECK(
      database_.sequences_ = (char **) malloc(
          sizeof(char *) * database_.num_of_sequences_));
  for (usint32_t i = 0; i < database_.num_of_sequences_; i++) {
    FREAD_CHECK(fread(&l, sizeof(usint32_t), 1, fin), 1);
    FREAD_CHECK(fread(proName, sizeof(char), l, fin), l);
    proName[l] = 0;
    database_.sequences_names_[i] = proName;
    FREAD_CHECK(fread(&l, sizeof(usint32_t), 1, fin), 1);
    //proDB_.vProSize.push_back(l);
    MEMORY_ALLOCATE_CHECK(
        database_.sequences_[i] = (char *) malloc(sizeof(char) * (l + 1)));
    FREAD_CHECK(fread(&(database_.sequences_[i][0]), sizeof(char), l, fin), l);
    database_.sequences_[i][l] = 0;
  }

  // kmer_nearest_
  FREAD_CHECK(
      fread(&(kmer_nearest_.size_of_counter_), sizeof(usint32_t), 1, fin), 1);
  FREAD_CHECK(fread(&(kmer_nearest_.size_of_index_), sizeof(usint64_t), 1, fin),
              1);
  MEMORY_ALLOCATE_CHECK(
      kmer_nearest_.counter_ = (usint64_t *) malloc(
          sizeof(usint64_t) * (kmer_nearest_.size_of_counter_ + 1)));
  MEMORY_ALLOCATE_CHECK(
      kmer_nearest_.index_ = (usint32_t *) malloc(
          sizeof(usint32_t) * kmer_nearest_.size_of_index_));
//  MEMORY_ALLOCATE_CHECK(
//      kmer_nearest_.score_ = (int * ) malloc(
//          sizeof(int) * kmer_nearest_.size_of_index_));
  FREAD_CHECK(
      fread(kmer_nearest_.counter_, sizeof(usint64_t),
            kmer_nearest_.size_of_counter_ + 1, fin),
      kmer_nearest_.size_of_counter_ + 1);
  FREAD_CHECK(
      fread(kmer_nearest_.index_, sizeof(usint32_t),
            kmer_nearest_.size_of_index_, fin),
      kmer_nearest_.size_of_index_);
//  FREAD_CHECK(
//      fread(kmer_nearest_.score_, sizeof(int), kmer_nearest_.size_of_index_,
//            fin),
//      kmer_nearest_.size_of_index_);

// kmer_position_
  FREAD_CHECK(fread(&(kmer_pos_.size_of_counter_), sizeof(usint32_t), 1, fin),
              1);
  FREAD_CHECK(fread(&(kmer_pos_.size_of_index_), sizeof(usint64_t), 1, fin), 1);
  MEMORY_ALLOCATE_CHECK(
      kmer_pos_.counter_ = (usint64_t *) malloc(
          sizeof(usint64_t) * (kmer_pos_.size_of_counter_ + 1)));
  MEMORY_ALLOCATE_CHECK(
      kmer_pos_.index_pro_id_ = (usint32_t *) malloc(
          sizeof(usint32_t) * kmer_pos_.size_of_index_));
  MEMORY_ALLOCATE_CHECK(
      kmer_pos_.index_pro_pos_ = (usint32_t *) malloc(
          sizeof(usint32_t) * kmer_pos_.size_of_index_));
  FREAD_CHECK(
      fread(kmer_pos_.counter_, sizeof(usint64_t),
            kmer_pos_.size_of_counter_ + 1, fin),
      kmer_pos_.size_of_counter_ + 1);
  FREAD_CHECK(
      fread(kmer_pos_.index_pro_id_, sizeof(usint32_t),
            kmer_pos_.size_of_index_, fin),
      kmer_pos_.size_of_index_);
  FREAD_CHECK(
      fread(kmer_pos_.index_pro_pos_, sizeof(usint32_t),
            kmer_pos_.size_of_index_, fin),
      kmer_pos_.size_of_index_);

  fclose(fin);
  //TIME_INFO(TestIndex(), "Write Test Index");
}
