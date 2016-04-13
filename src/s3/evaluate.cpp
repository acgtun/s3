/*
 *    This is the main function to evaluate different protein search tools.
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

/*
  ./evaluate -local test_10000_gap.local_alignment_all.out \
             -ssearch ssearch.out -fasta fasta36.out \
             -blastp blastp.out \
             -rapsearch2 rapsearch2.out.m8 \
             -ghostx ghostx.out \
             -diamond diamond.out \
             -sort \
             -K 1 \
             -M 10
 * */

#include "smithlab_os.hpp"
#include "OptionParser.hpp"

#include "evaluation.hpp"

string GetTimePre(const string& file_name) {
  time_t rawtime;
  struct tm * timeinfo;

  time(&rawtime);
  timeinfo = localtime(&rawtime);

  char pre[100] = { '\0' };
  sprintf(pre, "%04d-%02d-%02d-%02d-%02d-%02d_", 1900 + timeinfo->tm_year,
          1 + timeinfo->tm_mon, timeinfo->tm_mday, timeinfo->tm_hour,
          timeinfo->tm_min, timeinfo->tm_sec);

  string ret = pre;
  ret += file_name;
  return ret;
}


int main(int argc, const char* argv[]) {
  try {
    string localalign_file;
    string ssearch_file;
    string fasta_file;
    string blastp_file;
    string rapsearch2_file;
    string ghostx_file;
    string s3_file;
    string diamond_file;

    double evalue_threshold = 10;
    bool is_sort = true;
    int K = 1;
    int M = 10;
    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]),
                           "evaluate different protein sequence search tools",
                           "");
    opt_parse.add_opt("local", 'l', "local alignment results as ground truth",
                      true, localalign_file);
    opt_parse.add_opt("ssearch", 's', "ssearch", true, ssearch_file);
    opt_parse.add_opt("fasta", 'f', "fasta_file", true, fasta_file);
    opt_parse.add_opt("blastp", 'b', "blastp", true, blastp_file);
    opt_parse.add_opt("rapsearch2", 'r', "rapsearch2", true, rapsearch2_file);
    opt_parse.add_opt("ghostx", 'g', "ghostx", true, ghostx_file);
    opt_parse.add_opt("s3", 's', "S3", true, s3_file);
    opt_parse.add_opt("diamond", 'd', "diamond", true, diamond_file);

    opt_parse.add_opt("evalue", 'e', "evalue", false, evalue_threshold);
    opt_parse.add_opt("sort", 'p', "sort results", false, is_sort);
    opt_parse.add_opt("K", 'k', "top K", false, K);
    opt_parse.add_opt("M", 'm', "in top M", false, M);

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

    Evaluation localalign(localalign_file, NULL, evalue_threshold, is_sort, K,
                          M);
    Evaluation ssearch(ssearch_file, &localalign, evalue_threshold, is_sort, K,
                       M);
    Evaluation fasta(fasta_file, &localalign, evalue_threshold, is_sort, K, M);
    Evaluation blastp(blastp_file, &localalign, evalue_threshold, is_sort, K,
                      M);
    Evaluation rapsearch2(rapsearch2_file, &localalign, evalue_threshold,
                          is_sort, K, M);
    Evaluation ghostx(ghostx_file, &localalign, evalue_threshold, is_sort, K,
                      M);
    Evaluation s3(s3_file, &localalign, evalue_threshold, is_sort, K, M);
    Evaluation diamond(diamond_file, &localalign, evalue_threshold, is_sort, K,
                       M);

    ofstream fout(GetTimePre("accuracy.xlsx").c_str());
    cout << GetTimePre("accuracy.xlsx") << endl;

    for (double evalue = 1e-64; evalue <= 100; evalue *= 100) {
      fout << "\t" << evalue;
    }
    fout << endl;
    if (localalign_file.size() != 0) {
      localalign.ReadResults();
    }

    if (ssearch_file.size() != 0) {
      ssearch.ReadResults();
      fout << "ssearch36\t";
      ssearch.CompareTopKHitProtein(fout);
    }

    if (fasta_file.size() != 0) {
      fasta.ReadResults();
      fout << "fasta36\t";
      fasta.CompareTopKHitProtein(fout);
    }

    if (blastp_file.size() != 0) {
      blastp.ReadResults();
      fout << "blastp\t";
      blastp.CompareTopKHitProtein(fout);
    }

    if (rapsearch2_file.size() != 0) {
      rapsearch2.ReadResults();
      fout << "rapsearch2\t";
      rapsearch2.CompareTopKHitProtein(fout);
    }

    if (ghostx_file.size() != 0) {
      ghostx.ReadResults();
      fout << "ghostx\t";
      ghostx.CompareTopKHitProtein(fout);
    }

    if (s3_file.size() != 0) {
      s3.ReadResults();
      fout << "s3\t";
      s3.CompareTopKHitProtein(fout);
      //s3.CompareTopKHitProtein_INFO();
    }

    if (diamond_file.size() != 0) {
      diamond.ReadResults();
      fout << "diamond\t";
      diamond.CompareTopKHitProtein(fout);
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
