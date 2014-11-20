#include "./util/option.h"
#include "./index/index.h"
#include "./seqalign/local_alignment.h"
#include "query_search.h"

using local_alignment::LocalAlignment;

void Run(const FastaFile* queries, Evalue* evalue, QuerySearch* query_search) {
  cout << "Num of queries " << queries->num_of_sequences_ << endl;
  for (usint32_t i = 0; i < queries->num_of_sequences_; i++) {
    cout << "query " << i << endl;
    evalue->UpdateValues(strlen(queries->sequences_[i]));
    query_search->Search(queries->sequences_[i],
                         queries->sequences_names_[i].c_str());
  }
}

int main(int argc, const char* argv[]) {
  InitProgram(argc, argv);
  Index index;
  TIME_INFO(index.ReadIndex(), "Read Index");

  string query_file;
  Option::GetOption("-q", query_file);
  FastaFile queries(query_file);

  Evalue evalue(index.database_.num_of_characters_,
                index.database_.num_of_sequences_);

  QuerySearch query_search(&index, &evalue, queries.max_sequence_length_);
  TIME_INFO(Run(&queries, &evalue, &query_search), "Search");

  return 0;
}
