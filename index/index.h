#ifndef BUILDINDEX_H_
#define BUILDINDEX_H_

#include "./../util/bio_util.h"
#include "./../util/option.h"
#include "./../util/fasta_file.h"

#include "k_mer_position.h"
#include "k_mer_nearest.h"

#include <map>
#include <set>
#include <vector>
#include <string>

class Index {
 public:
  Index()
      : kmer_pos_(&database_),
        kmer_nearest_(&kmer_pos_) {
    string database_file;
    Option::GetOption("-db", database_file, "");
    Option::GetOption("-index", index_file_, "");
    if (Option::GetCommand() == "build") {
      if(database_file.size() == 0) {
        printf("Please use the parameter -db to set the database.");
        exit (EXIT_FAILURE);
      }
      index_file_ = Option::GetFileName(database_file);
      index_file_ += ".dbindex";
    } else if (Option::GetCommand() == "aligner") {
      if (index_file_.size() == 0) {
        printf("Please use the parameter -index to set the database.");
        exit (EXIT_FAILURE);
      }
      if (Option::GetFileExtension(index_file_) == "fasta") {
        ERROR_INFO(
            "Please use 'build' command to build index for the database first.");
        exit (EXIT_FAILURE);
      }
      if (Option::GetFileExtension(index_file_) != ".dbindex") {
        ERROR_INFO("The index file should be *.dbindex");
        exit (EXIT_FAILURE);
      }
    }

    Option::GetOption("-e", expected_hits_for_query_, 100);
    if (database_file.size() != 0) {
      database_.SetFilePath(database_file);
    }
  }

  ~Index() {
  }

  void BuildIndex();
  void ReadIndex();

 private:

  void WriteIndex();
  void GetThresholdScore();
  void TestIndex();

 public:
  FastaFile database_;
  KmerPosition kmer_pos_;
  KmerNearest kmer_nearest_;

  string index_file_;
  int expected_hits_for_query_;

}
;

#endif /* BUILDINDEX_H_ */
