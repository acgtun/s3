/*
 *    This is the main function to build index for protein database.
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

int main(int argc, const char* argv[]) {
  try {
    string database_file;
    string outfile;
    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]),
                           "build index for protein database", "");
    opt_parse.add_opt("database", 'd', "protein database file", true,
                      database_file);
    opt_parse.add_opt(
        "output", 'o',
        "output file name (the suffix of the file should be '.dbindex')", true,
        outfile);

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
    if (!is_valid_filename(outfile, "dbindex")) {
      fprintf(stderr, "The suffix of the output file should be '.dbindex'\n");
      return EXIT_FAILURE;
    }
    if (outfile.size() > 1000) {
      fprintf(stderr, "The output file name is too long, "
              "please select a shorter name\n");
      return EXIT_FAILURE;;
    }
    if (!leftover_args.empty()) {
      fprintf(stderr, "%s\n", opt_parse.help_message().c_str());
      return EXIT_SUCCESS;
    }
    /****************** END COMMAND LINE OPTIONS *****************/

    CProteinDB proteindb;
    HashTable hash_table;
    //////////////////////////////////////////////////////////////
    // READ DATABASE
    TIME_INFO(BuildProteinDB(database_file, proteindb), "READ DATABASE");

    // BUILD INDEX
    TIME_INFO(BuildKmerLocation(proteindb, hash_table), "BUILD KMER LOCATIONS");

    // WRITE IDNEX
    TIME_INFO(WriteIndex(proteindb, hash_table, outfile), "WRITE INDEX");
  } catch (const SMITHLABException &e) {
    fprintf(stderr, "%s\n", e.what().c_str());
    return EXIT_FAILURE;
  } catch (std::bad_alloc &ba) {
    fprintf(stderr, "ERROR: could not allocate memory\n");
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
