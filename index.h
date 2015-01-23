#ifndef INDEX_H_
#define INDEX_H_

#include "sdk.h"

#include <set>
#include <string>
#include <vector>
#include <utility>
#include <tr1/unordered_map>
#include <tr1/unordered_set>

using std::set;
using std::vector;
using std::string;
using std::make_pair;
using std::tr1::unordered_map;
using std::tr1::unordered_set;

struct DBLocation {
  DBLocation(const uint32_t& _protein_id, const uint32_t& _protein_pos)
      : protein_id(_protein_id),
        protein_pos(_protein_pos) {
  }
  DBLocation() {
    protein_id = 0;
    protein_pos = 0;
  }

  uint32_t protein_id;
  uint32_t protein_pos;
};

typedef std::tr1::unordered_map<uint32_t, vector<DBLocation> > KMERDBLOCATIONS;
typedef std::tr1::unordered_map<uint32_t, vector<uint32_t> > KMERNEIGHBORS;
typedef std::tr1::unordered_set<uint32_t> ITEM_SET;
typedef std::tr1::unordered_map<uint32_t, uint32_t> ITEM_COUNTING;
typedef std::tr1::unordered_map<uint32_t, set<uint32_t> > UNIQUE_ITEM_COUNTING;

struct CProtein {
  string name;
  uint32_t length;
  vector<char> sequence;
};

struct CProteinDB {
  vector<CProtein> proteins;
  uint64_t db_length;
  uint32_t num_of_proteins;
  uint32_t max_protein_length;
  uint32_t min_protein_length;
};

void ReadFASTAFile(const string& fasta_file, vector<string>& names,
                   vector<string>& seqs);
void BuildProteinDB(const string& database_file, CProteinDB& proteindb);
void BuildKmerLocation(const CProteinDB& proteindb,
                       KMERDBLOCATIONS& kmer_dblocations,
                       ITEM_SET& kmer_db_exist);
void BuildKmerNeighbors(KMERNEIGHBORS& kmer_neighbors,
                        const ITEM_SET& kmer_db_exist);

void WriteIndex(const CProteinDB& proteindb,
                const KMERDBLOCATIONS& kmer_dblocations,
                const KMERNEIGHBORS& kmer_neighbors);
void ReadIndex(CProteinDB& proteindb, KMERDBLOCATIONS& kmer_dblocations,
               KMERNEIGHBORS& kmer_neighbors);

#endif /* INDEX_H_ */
