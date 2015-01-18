#ifndef QUERY_SEARCH_H_
#define QUERYS_EARCH_H_

void GetTopProteinIDS(const KMERDBLOCATIONS& kmer_dblocations,
                      const KMERNEIGHBORS& kmer_neighbors,
                      const uint32_t& num_of_top_proteins,
                      const char* query_seq, ITEM_SET& protein_candidates_id);

#endif /* QUERY_SEARCH_H_ */
