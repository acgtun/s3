#include "sdk.h"
#include "option.h"
#include "index.h"

int main(int argc, const char* argv[]) {
  InitProgram(argc, argv);

  string database_file;
  Option::GetOption("-d", database_file);

  CProteinDB proteindb;
  KMERDBLOCATIONS kmer_dblocations;
  KMERNEIGHBORS kmer_neighbors;
  ITEM_SET kmer_db_exist;

  TIME_INFO(BuildProteinDB(database_file, proteindb), "READ DATABASE");
  TIME_INFO(BuildKmerLocation(proteindb, kmer_dblocations, kmer_db_exist),
            "BUILD KMER LOCATIONS");
  TIME_INFO(BuildKmerNeighbors(kmer_neighbors, kmer_db_exist),
            "BUILD KMER NEIGHBORS");
  TIME_INFO(WriteIndex(proteindb, kmer_dblocations, kmer_neighbors),
            "WRITE INDEX");

  return 0;
}
