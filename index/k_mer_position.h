#ifndef KMERPOSITION_H_
#define KMERPOSITION_H_

#include "./../util/fasta_file.h"
#include "./../util/option.h"

class KmerPosition {
 public:
  KmerPosition(const FastaFile* database);
  ~KmerPosition();

  void BuildKmerPosition();
  void CountKmerPosition();

 public:
  usint64_t* counter_;
  usint32_t* index_pro_id_;
  usint32_t* index_pro_pos_;
  usint32_t size_of_counter_;
  usint64_t size_of_index_;

 private:
  const FastaFile* database_;
  usint32_t MAX_BUCKET_SIZE;
};

#endif /* KMERPOSITION_H_ */
