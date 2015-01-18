#include "option.h"
#include "index.h"
#include "query_search.h"

int main(int argc, const char* argv[]) {
  InitProgram(argc, argv);

  CProteinDB proteindb;
  KMERDBLOCATIONS kmer_dblocations;
  KMERNEIGHBORS kmer_neighbors;

  TIME_INFO(ReadIndex(proteindb, kmer_dblocations, kmer_neighbors),
            "READ INDEX");

  string query_file;
  Option::GetOption("-q", query_file);

  vector<string> query_seqs, query_names;
  ReadFASTAFile(query_file, query_names, query_seqs);

  uint32_t num_of_top_proteins;
  Option::GetOption("-t", num_of_top_proteins, 50);

  uint32_t num_of_queries = query_seqs.size();
  for (uint32_t i = 0; i < num_of_queries; ++i) {
    ITEM_SET protein_candidates_id;
    GetTopProteinIDS(kmer_dblocations, kmer_neighbors, num_of_top_proteins,
                     query_seqs[i].c_str(), protein_candidates_id);
  }

  return 0;
}
