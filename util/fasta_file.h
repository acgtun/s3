#ifndef FASTAFILE_H_
#define FASTAFILE_H_

#include "bio_util.h"

class FastaFile {
 public:
  FastaFile(const string& file_path);
  void SetFilePath(const string& file_path);
  string GetFilePath() const;
  FastaFile();
  ~FastaFile();

 private:
  usint64_t GetLineFromString(const char* strVal, char* strRet);
  void GetNumOfSequences();
  void AnalyzeFile();
  void ReadFastaFile();

  string file_path_;  // Absolute path for the fasta file
  char* file_string_;  // To store characters of the whole file
  usint64_t file_size_;

 public:
  vector<string> sequences_names_;
  char** sequences_;
  usint32_t num_of_sequences_;
  usint64_t num_of_characters_;
  usint32_t max_sequence_length_;
};

#endif /* FASTAFILE_H_ */
