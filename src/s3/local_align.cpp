/*
 *    This is the main function for local alignment.
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

#include "smithlab_os.hpp"
#include "OptionParser.hpp"

#include "local_alignment.hpp"

#include "util.hpp"
#include "evalue.hpp"
#include "reference.hpp"

using local_alignment::LocalAlignment;

int main(int argc, const char* argv[]) {
  try {
    string database_file, query_file;
    double evalue_threshold = 10;
    int outfmt = 7;
    int gapopen = -11;
    int gapextension = -1;
    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]),
                           "local alignment using Smith-Waterman Algorithm",
                           "");
    opt_parse.add_opt("database", 'd', "protein or DNA database file",
                      true, database_file);
    opt_parse.add_opt("evalue", 'e', "evalue threshold", false, evalue_threshold);
    opt_parse.add_opt("format", 'f', "outout format", false, outfmt);

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
    vector<string> protein_names, protein_seqs;
    vector<string> query_names, query_seqs;
    ReadFASTAFile(database_file, protein_names, protein_seqs);
    ReadFASTAFile(query_file, query_names, query_seqs);

    uint64_t protein_dabase_length = 0, max_protein_length = 0;
    for (uint32_t i = 0; i < protein_seqs.size(); ++i) {
      protein_dabase_length += protein_seqs.size();
      max_protein_length =
          max_protein_length > protein_seqs[i].size() ?
              max_protein_length : protein_seqs[i].size();
    }

    uint32_t max_query_length = 0;
    for (uint32_t i = 0; i < query_seqs.size(); ++i) {
      max_query_length =
          max_query_length > query_seqs[i].size() ?
              max_query_length : query_seqs[i].size();
    }

    Evalue evalue(protein_dabase_length, protein_names.size(),
                  evalue_threshold);

    LocalAlignment local_align(&evalue, max_query_length, max_protein_length,
                               gapopen, gapextension);
    M8Results res;
    ofstream fout(string(query_file + "full_local_alignment.txt").c_str());
    for (uint32_t i = 0; i < query_seqs.size(); ++i) {
      cout << "query " << i << " " << query_names[i] << endl;
      evalue.UpdateValues(static_cast<uint32_t>(query_seqs[i].size()));
      vector<M8Results> aligned_results;
      for (uint32_t j = 0; j < protein_seqs.size(); ++j) {
        if (local_align.RunLocalAlignment(query_seqs[i], protein_seqs[j],
                                          res)) {
          cout << "protein: " << j << " " << protein_names[j] << endl;
          res.protein_name = protein_names[j];
          aligned_results.push_back(res);
        }
      }
      DisplayResults(query_names[i], database_file, aligned_results,
                     aligned_results.size(), outfmt, fout);
    }
    fout.close();
  } catch (const SMITHLABException &e) {
    fprintf(stderr, "%s\n", e.what().c_str());
    return EXIT_FAILURE;
  } catch (std::bad_alloc &ba) {
    fprintf(stderr, "ERROR: could not allocate memory\n");
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
