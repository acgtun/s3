/*
 *    This is the main function of S3.
 *
 *    Copyright (C) 2016 University of Southern California
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

#include "option.h"
#include "index.h"
#include "query_search.h"
#include "local_alignment.h"

using local_alignment::LocalAlignment;

void Matching(const CProteinDB& proteindb, const HashTable& hash_table,
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
  Option::GetOption("-t", num_of_top_proteins, 100);

  uint32_t outfmt;
  Option::GetOption("-f", outfmt, 8);

  string output_file;
  Option::GetOption("-o", output_file);
  ofstream fout(output_file.c_str());

  uint32_t num_of_protein_candidates = 0;
  vector < uint32_t > protein_candidates_id(num_of_top_proteins);

  Evalue evalue(proteindb.db_length, proteindb.num_of_proteins,
                evalue_threshold);
  LocalAlignment local_alignment(&evalue, max_query_length,
                                 proteindb.max_protein_length);
  INFO("START MAPPING...");
  uint32_t num_of_queries = query_seqs.size();
  for (uint32_t i = 0; i < num_of_queries; ++i) {
    if (query_seqs[i].size() < HASHLEN)
      continue;
    GetTopProteinIDS(hash_table, query_seqs[i].c_str(), num_of_top_proteins,
                     protein_candidates_id, num_of_protein_candidates);
    evalue.UpdateValues(query_seqs[i].size());
    M8Results res;
    vector < M8Results > aligned_results;
    for (uint32_t j = 0; j < num_of_protein_candidates; ++j) {
      if (local_alignment.RunLocalAlignment(
          query_seqs[i], proteindb.proteins[protein_candidates_id[j]].sequence,
          res)) {
        res.protein_name = proteindb.proteins[protein_candidates_id[j]].name;
        aligned_results.push_back(res);
      }
    }
    DisplayResults(query_names[i], "", aligned_results, aligned_results.size(),
                   outfmt, fout);
  }
  fout.close();
}

int main(int argc, const char* argv[]) {
  InitProgram(argc, argv);

  CProteinDB proteindb;
  HashTable hash_table;

  TIME_INFO(ReadIndex(proteindb, hash_table), "READ INDEX");

  string query_file;
  Option::GetOption("-q", query_file);

  vector<string> query_seqs, query_names;
  ReadFASTAFile(query_file, query_names, query_seqs);
  INFO("THE NUMBER OF QUERIES IS", query_seqs.size());

  TIME_INFO(Matching(proteindb, hash_table, query_seqs, query_names),
            "PROTEIN MATCHING");

  return 0;
}
