#include "evaluation.h"
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, const char* argv[]) {
  InitProgram(argc, argv);
  string localalign_file, ssearch_file, fasta_file;
  string blastp_file, rapsearch2_file, ghostx_file, s3_file, s3_local_file;

  Option::GetOption("-local", localalign_file, "");
  Option::GetOption("-ssearch", ssearch_file, "");
  Option::GetOption("-fasta", fasta_file, "");
  Option::GetOption("-blastp", blastp_file, "");
  Option::GetOption("-rapsearch2", rapsearch2_file, "");
  Option::GetOption("-ghostx", ghostx_file, "");
  Option::GetOption("-s3", s3_file, "");
  Option::GetOption("-s3local", s3_local_file, "");

  Evaluation localalign(localalign_file, NULL);
  Evaluation ssearch(ssearch_file, &localalign);
  Evaluation fasta(fasta_file, &localalign);
  Evaluation blastp(blastp_file, &localalign);
  Evaluation rapsearch2(rapsearch2_file, &localalign);
  Evaluation ghostx(ghostx_file, &localalign);
  Evaluation s3(s3_file, &localalign);
  Evaluation s3_local(s3_local_file, &localalign);

  ofstream fout(Option::AddPreToFileName("accuracy.txt").c_str());
  cout << Option::AddPreToFileName("accuracy.txt") << endl;

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
    s3.CompareTopKHitProtein_INFO();
  }

  if (s3_local_file.size() != 0) {
    s3_local.ReadResults();
    fout << "s3_local\t";
    s3_local.CompareTopKHitProtein(fout);
  }

  fout.close();

  return 0;
}
