#include "k_mer_position.h"

KmerPosition::KmerPosition(const FastaFile* database)
    : database_(database) {
  counter_ = NULL;
  index_pro_id_ = NULL;
  index_pro_pos_ = NULL;
  size_of_counter_ = 0;
  size_of_index_ = 0;
  Option::GetOption("-bucket", MAX_BUCKET_SIZE, 5000);
}

KmerPosition::~KmerPosition() {
  free(counter_);
  free(index_pro_id_);
  free(index_pro_pos_);
}

void KmerPosition::CountKmerPosition() {
  INFO("Count Kmer Position...");
  size_of_counter_ = GetPower(20, HASHAALEN);
  INFO("The size of Kmer Position Counter is", size_of_counter_);
  MEMORY_ALLOCATE_CHECK(
      counter_ = (usint64_t * ) malloc(
          sizeof(usint64_t) * (size_of_counter_ + 1)));
  memset(counter_, 0x00, sizeof(usint64_t) * (size_of_counter_ + 1));

  usint32_t size = 0, hashValue = 0;
  for (usint32_t i = 0; i < database_->num_of_sequences_; i++) {
    //cout << "i = " << i << " " << database_->num_of_sequences_ << endl;
    size = strlen(database_->sequences_[i]);
    if (size < HASHAALEN) {
      INFO("The length of the protein is less than HASHAALEN, ignore it.",
           database_->sequences_names_[i].c_str());
      continue;
    }
    size = size - HASHAALEN;
    for (usint32_t j = 0; j <= size; j++) {
      hashValue = GetHashValue(&(database_->sequences_[i][j]));
      counter_[hashValue]++;
    }
  }

  for (usint32_t i = 1; i <= size_of_counter_; i++) {
    if (counter_[i] > MAX_BUCKET_SIZE) {
      counter_[i] = 0;
    }
  }

  for (usint32_t i = 1; i <= size_of_counter_; i++) {
    counter_[i] += counter_[i - 1];
  }
  size_of_index_ = counter_[size_of_counter_];
  INFO("The size of Kmer Position Index is", size_of_index_);

  for (usint32_t i = size_of_counter_ - 1; i >= 1; i--) {
    counter_[i] = counter_[i - 1];
  }
  counter_[0] = 0;
}

void KmerPosition::BuildKmerPosition() {
  INFO("Build Kmer Position...");
  MEMORY_ALLOCATE_CHECK(
      index_pro_id_ = (usint32_t * ) malloc(
          sizeof(usint32_t) * size_of_index_));
  MEMORY_ALLOCATE_CHECK(
      index_pro_pos_ = (usint32_t * ) malloc(
          sizeof(usint32_t) * size_of_index_));
  usint32_t size = 0, hashValue = 0;
  for (usint32_t i = 0; i < database_->num_of_sequences_; i++) {
    size = strlen(database_->sequences_[i]);
    if (size < HASHAALEN) {
      INFO("The length of the protein is less than 6, ignore it.",
           database_->sequences_names_[i].c_str());
      continue;
    }

    size = size - HASHAALEN;
    for (usint32_t j = 0; j <= size; j++) {
      hashValue = GetHashValue(&(database_->sequences_[i][j]));
      if (counter_[hashValue] == counter_[hashValue + 1])
        continue;
      index_pro_id_[counter_[hashValue]] = i;
      index_pro_pos_[counter_[hashValue]] = j;
      counter_[hashValue]++;
    }
  }

  for (usint32_t i = size_of_counter_ - 1; i >= 1; i--) {
    counter_[i] = counter_[i - 1];
  }
  counter_[0] = 0;
}
