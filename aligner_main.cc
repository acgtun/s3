#include "option.h"
#include "index.h"

int main(int argc, const char* argv[]) {
  InitProgram(argc, argv);

  CProteinDB proteindb;
  KMERDBLOCATIONS kmer_dblocations;
  KMERNEIGHBORS kmer_neighbors;

  TIME_INFO(ReadIndex(proteindb, kmer_dblocations, kmer_neighbors),
            "READ INDEX");

  return 0;
}
