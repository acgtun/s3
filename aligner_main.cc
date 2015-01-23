#include "option.h"
#include "index.h"
#include "query_search.h"
#include "local_alignment.h"

using local_alignment::LocalAlignment;

void Matching(const CProteinDB& proteindb,
              const KMERDBLOCATIONS& kmer_dblocations,
              const KMERNEIGHBORS& kmer_neighbors,
              const vector<string>& query_seqs,
              const vector<string>& query_names) {
  uint32_t max_query_length = 0;
  for (uint32_t i = 0; i < query_seqs.size(); ++i) {
    max_query_length =
        max_query_length > query_seqs[i].size() ?
            max_query_length : query_seqs[i].size();
  }
  INFO("THE MAXIMAL LENGTH OF A QUERY IS", max_query_length);

  uint32_t num_of_top_proteins;
  Option::GetOption("-t", num_of_top_proteins, 50);

  uint32_t outfmt;
  Option::GetOption("-f", outfmt, 8);

  string output_file;
  Option::GetOption("-o", output_file);
  ofstream fout(output_file.c_str());

  vector<uint32_t> protein_count(proteindb.num_of_proteins);

  uint32_t num_of_protein_candidates = 0;
  vector<uint32_t> protein_candidates_id(num_of_top_proteins);

  Evalue evalue(proteindb.db_length, proteindb.num_of_proteins);
  LocalAlignment local_alignment(&evalue, max_query_length,
                                 proteindb.max_protein_length);

  INFO("START MAPPING...");
  vector<M8Results> aligned_results(num_of_top_proteins);
  uint32_t num_of_queries = query_seqs.size();
  for (uint32_t i = 0; i < num_of_queries; ++i) {
    if(query_seqs[i].size() < HASHLEN) continue;
   // cout << "read " << i << " ";
    GetTopProteinIDS(proteindb.num_of_proteins, kmer_dblocations,
                     kmer_neighbors, num_of_top_proteins, protein_count,
                     query_seqs[i].c_str(), protein_candidates_id,
                     num_of_protein_candidates);
    evalue.UpdateValues(query_seqs[i].size());
    M8Results res;
    uint32_t num_of_results = 0;
    for (uint32_t j = 0; j < num_of_protein_candidates; ++j) {
      //cout << i << " " << protein_candidates_id[j] << endl;
      //cout << query_seqs[i] << endl;
      //for (uint32_t k = 0;
         // k < proteindb.proteins[protein_candidates_id[j]].sequence.size();
          //++k) {
        //cout << proteindb.proteins[protein_candidates_id[j]].sequence[k];
     // }
     // cout << endl;
      //cout << "read " << i << endl;
      if (local_alignment.RunLocalAlignment(
          query_seqs[i], proteindb.proteins[protein_candidates_id[j]].sequence,
          res)) {
        res.protein_id = protein_candidates_id[j];
        aligned_results[num_of_results++] = res;
      }
    }
    DisplayResults(proteindb, query_names[i], "", aligned_results,
                   num_of_results, outfmt, fout);
  }
  fout.close();
}

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
  INFO("THE NUMBER OF QUERIES IS", query_seqs.size());

  TIME_INFO(
      Matching(proteindb, kmer_dblocations, kmer_neighbors, query_seqs,
               query_names),
      "PROTEIN MATCHING");

  return 0;
}
