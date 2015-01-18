#include "index.h"
#include "option.h"
#include "bio_util.h"
#include "query_search.h"

#include <queue>
#include <functional>
using std::greater;
using std::priority_queue;

void TopProteinsThreshold(const ITEM_COUNTING& protein_matched_kmer_count,
                          const uint32_t& num_of_top_proteins,
                          uint32_t& protein_count_threshold) {
  if (protein_matched_kmer_count.size() <= num_of_top_proteins) {
    protein_count_threshold = 0;
    return;
  }
  priority_queue<uint32_t, vector<uint32_t>, greater<uint32_t> > top_proteins;
  for (ITEM_COUNTING::const_iterator it = protein_matched_kmer_count.begin();
      it != protein_matched_kmer_count.end(); ++it) {
    if (top_proteins.size() < num_of_top_proteins) {
      top_proteins.push(it->second);
    } else if (top_proteins.size() == num_of_top_proteins
        && it->second > top_proteins.top()) {
      top_proteins.pop();
      top_proteins.push(it->second);
    }
  }
  protein_count_threshold = top_proteins.top();
}

void GetTopProteinIDS(const KMERDBLOCATIONS& kmer_dblocations,
                      const KMERNEIGHBORS& kmer_neighbors,
                      const uint32_t& num_of_top_proteins,
                      const char* query_seq, ITEM_SET& protein_candidates_id) {
  ITEM_COUNTING protein_matched_kmer_count;
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
        protein_matched_kmer_count[it_kmer->second[p].protein_id]++;
      }
    }
  }
  uint32_t protein_count_threshold = 0;
  TopProteinsThreshold(protein_matched_kmer_count, num_of_top_proteins,
                       protein_count_threshold);

  for (ITEM_COUNTING::const_iterator it = protein_matched_kmer_count.begin();
      it != protein_matched_kmer_count.end(); ++it) {
    if (it->second >= protein_count_threshold) {
      protein_candidates_id.insert(it->first);
    }
  }
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
