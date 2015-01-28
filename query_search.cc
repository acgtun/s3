#include "query_search.h"

#include <algorithm>

using std::make_heap;
using std::pop_heap;

void TopProteinsThreshold(vector<uint32_t>& protein_count,
                          const uint32_t& num_of_top_proteins,
                          uint32_t& protein_count_threshold) {
  make_heap(protein_count.begin(), protein_count.end());
  for (unsigned j = 0; j < num_of_top_proteins; j++) {
    pop_heap(protein_count.begin(), protein_count.end() - j);
  }

  protein_count_threshold = protein_count.front();
}

void GetTopProteinIDS(const uint32_t& num_of_proteins,
                      const KMERDBLOCATIONS& kmer_dblocations,
                      const KMERNEIGHBORS& kmer_neighbors,
                      const uint32_t& num_of_top_proteins,
                      vector<uint32_t>& protein_count, const char* query_seq,
                      vector<uint32_t>& protein_candidates_id,
                      uint32_t& num_of_protein_candidates) {
  memset(&(protein_count[0]), 0, sizeof(uint32_t) * num_of_proteins);
  //ITEM_COUNTING protein_matched_kmer_count;
  uint32_t num_of_seed = strlen(query_seq) - HASHLEN + 1;
  // cout << "xfsdax" << endl;
  for (uint32_t j = 0; j < num_of_seed; j++) {
    /* (1) Get Hash Value of the Kmer (here Kmer is seed) */
    uint32_t hash_value = Kmer2Integer(&(query_seq[j]));
    //cout << j << " " << num_of_seed << " " << hash_value << endl;
    /* (2) Get Nearest Kmers from the kmer_nearest table */
    KMERNEIGHBORS::const_iterator it = kmer_neighbors.find(hash_value);
    if (it == kmer_neighbors.end())
      continue;
    //cout << "it->second.size() = " << it->second.size() << endl;
    for (uint32_t c = 0; c < it->second.size(); ++c) {
      /* (3) Get the Kmer */
      uint32_t kmer = it->second[c];
      KMERDBLOCATIONS::const_iterator it_kmer = kmer_dblocations.find(kmer);
      if (it_kmer == kmer_dblocations.end())
        continue;

      /* (4) Get all the protein database position which occurs the Kmer */
      for (uint32_t p = 0; p < it_kmer->second.size(); ++p) {
        //cout << it_kmer->second[p].protein_id << endl;
        protein_count[it_kmer->second[p].protein_id]++;
      }
    }
  }
  // cout << "xx" << endl;
  uint32_t protein_count_threshold = 0;
  TopProteinsThreshold(protein_count, num_of_top_proteins,
                       protein_count_threshold);

  //uint32_t cnt = 0;
  num_of_protein_candidates = 0;
  for (uint32_t i = 0; i < protein_count.size(); ++i) {
    //if(protein_count[i] != 0) cnt++;
    if (protein_count[i] > protein_count_threshold) {
      protein_candidates_id[num_of_protein_candidates++] = i;
      //cout << num_of_protein_candidates << endl;
    }
  }
  //cerr << "cnt = " << cnt << endl;
}

void QuerySearch(const CProteinDB& proteindb,
                 const KMERDBLOCATIONS& kmer_dblocations,
                 const KMERNEIGHBORS& kmer_neighbors, const char* query_seq,
                 const ITEM_SET& protein_candidates_id) {
  KMERDBLOCATIONS protein_hits;
  uint32_t num_of_seed = strlen(query_seq) - HASHLEN + 1;
  for (uint32_t j = 0; j < num_of_seed; j++) {
    /* (1) Get Hash Value of the Kmer (here Kmer is seed) */
    uint32_t hash_value = Kmer2Integer(&(query_seq[j]));

    /* (2) Get Nearest Kmers from the kmer_nearest table */
    KMERNEIGHBORS::const_iterator it = kmer_neighbors.find(hash_value);
    if (it == kmer_neighbors.end())
      continue;

    for (uint32_t c = 0; c < it->second.size(); ++c) {
      /* (3) Get the Kmer */
      uint32_t kmer = it->second[c];
      KMERDBLOCATIONS::const_iterator it_kmer = kmer_dblocations.find(kmer);
      if (it_kmer == kmer_dblocations.end())
        continue;

      /* (4) Get all the protein database position which occurs the Kmer */
      for (uint32_t p = 0; p < it_kmer->second.size(); ++p) {
        if (protein_candidates_id.find(it_kmer->second[p].protein_id)
            != protein_candidates_id.end()) {
          protein_hits[it_kmer->second[p].protein_id].push_back(
              DBLocation(j, it_kmer->second[p].protein_pos));
        }
      }
    }
  }
}

