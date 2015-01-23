#ifndef QUERY_SEARCH_H_
#define QUERYS_EARCH_H_

#include "sdk.h"
#include "index.h"
#include "option.h"
#include "bio_util.h"

#include <vector>
#include <fstream>

using std::vector;
using std::ofstream;

void GetTopProteinIDS(const uint32_t& num_of_proteins,
                      const KMERDBLOCATIONS& kmer_dblocations,
                      const KMERNEIGHBORS& kmer_neighbors,
                      const uint32_t& num_of_top_proteins,
                      vector<uint32_t>& protein_count, const char* query_seq,
                      vector<uint32_t>& protein_candidates_id,
                      uint32_t& num_of_protein_candidates);

void DisplayResults(const CProteinDB& proteindb, const string& query_name,
                    const string& database_file,
                    const vector<M8Results>& aligned_results,
                    const uint32_t& num_of_results, const int& outfmt,
                    ofstream& fout);

#endif /* QUERY_SEARCH_H_ */
