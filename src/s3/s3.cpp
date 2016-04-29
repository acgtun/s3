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

#include "smithlab_os.hpp"
#include "OptionParser.hpp"

#include "reference.hpp"
#include "query_search.hpp"
#include "local_alignment.hpp"

using local_alignment::LocalAlignment;

void Matching(const CProteinDB& proteindb, const HashTable& hash_table,
              Evalue& evalue, LocalAlignment& local_alignment,
              const uint32_t& num_of_top_proteins, const int& outfmt,
              const string& output_file, const vector<string>& query_seqs,
              const vector<string>& query_names) {
  uint32_t num_of_candidates = 0;
  vector < uint32_t > candidates_id(num_of_top_proteins);
  fprintf(stderr, "START MAPPING...");
  ofstream fout(output_file.c_str());
  for (uint32_t i = 0; i < query_seqs.size(); ++i) {
    if (query_seqs[i].size() < HASHLEN)
      continue;
    GetTopProteinIDS(hash_table, query_seqs[i].c_str(), num_of_top_proteins,
                     candidates_id, num_of_candidates);
    evalue.UpdateValues(query_seqs[i].size());
    M8Results res;
    vector<M8Results> aligned_results;
    for (uint32_t j = 0; j < num_of_candidates; ++j) {
      if (local_alignment.RunLocalAlignment(
          query_seqs[i], proteindb.proteins[candidates_id[j]].sequence,
          res)) {
        res.protein_name = proteindb.proteins[candidates_id[j]].name;
        aligned_results.push_back(res);
      }
    }
    DisplayResults(query_names[i], "", aligned_results, aligned_results.size(),
                   outfmt, fout);
  }
  fout.close();
}

int main(int argc, const char* argv[]) {
  try {
    string index_file;
    string query_file;
    string output_file;

    double evalue_threshold = 10;

    // number of top proteins
    int top_k = 100;
    int outfmt = 8;
    int gapopen = -11;
    int gapextension = -1;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]),
                           "S3: protein sequence search tool", "");
    opt_parse.add_opt(
        "index",
        'i',
        "index file created by makedb command \
         (the suffix of the index file should be '.dbindex')",
        true, index_file);
    opt_parse.add_opt(
        "query",
        'q',
        "protein query sequences \
         (the suffix of read files should be '.fastq' or '.fq')",
        true, query_file);
    opt_parse.add_opt("output", 'o', "output file name", true, output_file);
    opt_parse.add_opt("evalue", 'e', "evalue threshold", false,
                      evalue_threshold);
    opt_parse.add_opt("format", 'f', "outout format", false, outfmt);
    opt_parse.add_opt("topk", 'k', "maximum number of results for each query",
                      false, top_k);

    vector<string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (argc == 1 || opt_parse.help_requested()) {
      fprintf(stderr, "%s\n", opt_parse.help_message().c_str());
      return EXIT_SUCCESS;
    }
    if (opt_parse.about_requested()) {
      fprintf(stderr, "%s\n", opt_parse.about_message().c_str());
      return EXIT_SUCCESS;
    }
    if (opt_parse.option_missing()) {
      fprintf(stderr, "%s\n", opt_parse.option_missing_message().c_str());
      return EXIT_SUCCESS;
    }
    if (!leftover_args.empty()) {
      fprintf(stderr, "%s\n", opt_parse.help_message().c_str());
      return EXIT_SUCCESS;
    }
    /****************** END COMMAND LINE OPTIONS *****************/

    CProteinDB proteindb;
    HashTable hash_table;
    TIME_INFO(ReadIndex(proteindb, hash_table, index_file), "READ INDEX");

    vector<string> query_seqs, query_names;
    ReadFASTAFile(query_file, query_names, query_seqs);
    uint32_t max_query_length = 0;
    for (uint32_t i = 0; i < query_seqs.size(); ++i) {
      max_query_length =
          max_query_length > query_seqs[i].size() ?
              max_query_length : query_seqs[i].size();
    }
    fprintf(stderr, "THE MAXIMAL LENGTH OF A QUERY IS %u\n.", max_query_length);

    Evalue evalue(proteindb.db_length, proteindb.num_of_proteins,
                  evalue_threshold);
    LocalAlignment local_alignment(&evalue, max_query_length,
                                   proteindb.max_protein_length, gapopen,
                                   gapextension);
    TIME_INFO(
        Matching(proteindb, hash_table, evalue, local_alignment, top_k, outfmt,
                 output_file, query_seqs, query_names),
        "PROTEIN MATCHING");
  } catch (const SMITHLABException &e) {
    fprintf(stderr, "%s\n", e.what().c_str());
    return EXIT_FAILURE;
  } catch (std::bad_alloc &ba) {
    fprintf(stderr, "ERROR: could not allocate memory\n");
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
