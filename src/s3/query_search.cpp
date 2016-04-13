/*
 *    This file contains functions for searching proteins.
 *
 *    Copyright (C) 2015 University of Southern California
 *
 *    Authors: Haifeng Chen and Ting Chen
 *
 *    This file is part of S3.
 *
 *    S3 is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    S3 is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with S3.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "query_search.hpp"
#include <queue>
#include <algorithm>

using std::priority_queue;

void TopProteinsThreshold(const unordered_map<uint32_t, uint32_t>& protein_hits,
                          const uint32_t& num_of_top_proteins,
                          uint32_t& protein_count_threshold) {
  priority_queue<uint32_t, vector<uint32_t>, std::greater<uint32_t> > top_protein;
  for (unordered_map<uint32_t, uint32_t>::const_iterator it =
      protein_hits.begin(); it != protein_hits.end(); ++it) {
    top_protein.push(it->second);
    if (top_protein.size() == num_of_top_proteins + 1) {
      top_protein.pop();
    }
  }

  if (top_protein.size() < num_of_top_proteins) {
    protein_count_threshold = 0;
  } else {
    protein_count_threshold = top_protein.top();
  }
}

void GetTopProteinIDS(const HashTable& hash_table,
                      const char* query_seq,
                      const uint32_t& num_of_top_proteins,
                      vector<uint32_t>& protein_candidates_id,
                      uint32_t& num_of_protein_candidates) {
  unordered_map<uint32_t, uint32_t> protein_hits;
  uint32_t num_of_seed = strlen(query_seq) - HASHLEN + 1;
  for (uint32_t j = 0; j < num_of_seed; j++) {
    uint32_t hash_value = Kmer2Integer(&(query_seq[j]));
    HashTable::const_iterator it = hash_table.find(hash_value);
    if (it == hash_table.end())
      continue;
    for (uint32_t p = 0; p < it->second.size(); ++p) {
      protein_hits[it->second[p].first]++;
    }
  }
  uint32_t protein_count_threshold = 0;
  TopProteinsThreshold(protein_hits, num_of_top_proteins,
                       protein_count_threshold);
  num_of_protein_candidates = 0;
  for (unordered_map<uint32_t, uint32_t>::iterator it = protein_hits.begin();
      it != protein_hits.end(); ++it) {
    if (it->second > protein_count_threshold) {
      protein_candidates_id[num_of_protein_candidates++] = it->first;
    }
  }
}

