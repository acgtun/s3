#include "option.h"
#include "index.h"
#include "bio_util.h"
#include "nearest_kmer.h"

#include <limits>
#include <string>
#include <fstream>

using std::string;
using std::ifstream;
using nearest_kmer::KNearestNeighbor;

void ReadDBFile(const string& database_file, vector<string>& protein_names,
                vector<string>& protein_seqs) {
  INFO("READ DATABASE FILE:", database_file.c_str());
  ifstream fin(database_file.c_str());
  if (!fin.good()) {
    ERROR_INFO("DATABASE FILE OPEN ERROR");
  }
  string sequence, line;
  while (getline(fin, line)) {
    if (line[0] == '>') {
      if (sequence.size() != 0) {
        protein_seqs.push_back(sequence);
        sequence.clear();
      }
      uint32_t space_pos = line.find_first_of(' ');
      if (space_pos == string::npos) {
        protein_names.push_back(line.substr(1));
      } else {
        protein_names.push_back(line.substr(1, space_pos - 1));
      }
      continue;
    }
    for (uint32_t i = 0; i < line.size(); ++i) {
      if (AA20.find_first_of(line[i]) != string::npos) {
        sequence += toupper(line[i]);
      } else if (isalpha(line[i])) {
        sequence += AA20[rand() % 20];  //todo
      }
    }
  }
  if (sequence.size() != 0) {
    protein_seqs.push_back(sequence);
  }
  fin.close();

  //////////////////////////////////////////////////////
  FILE * fout = fopen("proteins.txt", "w");
  for (uint32_t i = 0; i < protein_seqs.size(); i++) {
    fprintf(fout, "%d:%s\n", i, protein_names[i].c_str());
    fprintf(fout, "%s\n", protein_seqs[i].c_str());
  }
  fclose(fout);
  ///////////////////////////////////////////////////////
}

void BuildProteinDB(const string& database_file, CProteinDB& proteindb) {
  vector<string> protein_names, protein_seqs;
  ReadDBFile(database_file, protein_names, protein_seqs);

  proteindb.db_length = 0;
  proteindb.num_of_proteins = 0;
  proteindb.max_protein_length = 0;
  proteindb.min_protein_length = std::numeric_limits<unsigned int>::max();
  proteindb.proteins.resize(protein_seqs.size());
  for (uint32_t i = 0; i < protein_seqs.size(); ++i) {
    CProtein protein;
    protein.name = protein_names[i];
    protein.length = protein_seqs[i].size();
    protein.sequence.resize(protein.length);
    for (uint32_t j = 0; j < protein.length; ++j) {
      protein.sequence[j] = protein_seqs[i][j];
    }
    proteindb.proteins[i] = protein;
    proteindb.db_length += protein.length;
    proteindb.num_of_proteins++;
    proteindb.max_protein_length =
        proteindb.max_protein_length > protein.length ?
            proteindb.max_protein_length : protein.length;
    proteindb.min_protein_length =
        proteindb.min_protein_length < protein.length ?
            proteindb.min_protein_length : protein.length;
  }
}

void BuildKmerLocation(const CProteinDB& proteindb,
                       KMERDBLOCATIONS& kmer_dblocations,
                       KMER_DB_EXIST& kmer_db_exist) {
  INFO("BUILD KMER LOCATIONS...");
  for (uint32_t i = 0; i < proteindb.num_of_proteins; ++i) {
    const CProtein& protein = proteindb.proteins[i];
    if (protein.length < HASHLEN)
      continue;
    uint32_t size = protein.length - HASHLEN;
    for (uint32_t j = 0; j < size; ++j) {
      uint32_t hash_value = Kmer2Integer(&(protein.sequence[j]));
      kmer_dblocations[hash_value].push_back(DBLocation(i, j));
      kmer_db_exist.insert(hash_value);
    }
  }

#ifdef TEST
  for (KMERDBLOCATIONS::iterator it = kmer_dblocations.begin();
      it != kmer_dblocations.end(); ++it) {
    cout << it->first << endl;
    for (uint32_t i = 0; i < it->second.size(); ++i) {
      cout << " " << it->second[i].protein_id << " "
      << it->second[i].protein_pos << endl;
    }
  }
#endif
}

void BuildKmerNeighbors(KMERNEIGHBORS& kmer_neighbors,
                        const KMER_DB_EXIST& kmer_db_exist) {
  INFO("BUILD KMER NEIGHBORS...");
  uint32_t num_of_kmer = static_cast<uint32_t>(pow(ALPHABETSIZE, HASHLEN));
  uint32_t max_num_of_neighbors;
  Option::GetOption("-m", max_num_of_neighbors, 100);
  KNearestNeighbor k_nearest_neighbor;
  for (uint32_t i = 0; i < num_of_kmer; ++i) {
    //cout << i << " || ";
    k_nearest_neighbor.UpdateWeight(i);
    uint32_t path_id = 0, candidate = 0;
    double raw_score = 0.0;
    bool bpath_exist = true;
    while (bpath_exist && path_id < max_num_of_neighbors) {
      bpath_exist = k_nearest_neighbor.FindCandidate(candidate, path_id,
                                                     raw_score);
      if (bpath_exist && kmer_db_exist.find(candidate) != kmer_db_exist.end()) {
        kmer_neighbors[i].push_back(candidate);
        //cout << " " << candidate;
      }
      path_id++;
    }
    //cout << endl;
  }

#ifdef TEST
  for (KMERNEIGHBORS::iterator it = kmer_neighbors.begin();
      it != kmer_neighbors.end(); ++it) {
    cout << it->first << " || ";
    for (uint32_t i = 0; i < it->second.size(); ++i) {
      cout << " " << it->second[i];
    }
    cout << endl;
  }
#endif
}
