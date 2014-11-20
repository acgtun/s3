#include "./../util/option.h"
#include <string.h>
#include "./../util/evalue.h"
#include "./../util/fasta_file.h"
#include "local_alignment.h"
#include "./../util/bio_util.h"

using local_alignment::LocalAlignment;

int main(int argc, const char* argv[]) {
  InitProgram(argc, argv);
  string database_file, query_file;
  int outfmt;
  Option::GetOption("-db", database_file);
  Option::GetOption("-query_file", query_file);
  Option::GetOption("-outfmt", outfmt, 7);

  FastaFile database(database_file);
  FastaFile queries(query_file);
  Evalue evalue(database.num_of_characters_, database.num_of_sequences_);

  cout << queries.max_sequence_length_ << endl;
  cout << database.max_sequence_length_ << endl;

  LocalAlignment local_align(&evalue, queries.max_sequence_length_,
                             database.max_sequence_length_);
  M8Results res;
  ofstream fout(string(query_file + "full_local_alignment.txt").c_str());
  for (usint32_t i = 0; i < queries.num_of_sequences_; ++i) {
    cout << "query " << i << " " << queries.sequences_names_[i] << endl;
    evalue.UpdateValues(strlen(queries.sequences_[i]));
    vector<M8Results> aligned_results;
    for (usint32_t j = 0; j < database.num_of_sequences_; ++j) {
      if (local_align.RunLocalAlignment(queries.sequences_[i],
                                        database.sequences_[j],
                                        strlen(queries.sequences_[i]),
                                        strlen(database.sequences_[j]), res)) {
        cout << "protein: " << j << " " << database.sequences_names_[j] << endl;
        res.protein_name = database.sequences_names_[j];
        aligned_results.push_back(res);
      }
    }
    DisplayResults(queries.sequences_names_[i].c_str(), database_file.c_str(),
                   aligned_results, outfmt, fout);
  }
  fout.close();

  cout << "hree" << endl;
  return 0;
}
