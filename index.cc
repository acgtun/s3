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
  INFO("READ DATABASE FILE", database_file.c_str());
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
  /* Count Bucket Size */
  std::tr1::unordered_map<uint32_t, uint32_t> kmer_count;
  for (uint32_t i = 0; i < proteindb.num_of_proteins; ++i) {
    const CProtein& protein = proteindb.proteins[i];
    if (protein.length < HASHLEN)
      continue;
    uint32_t size = protein.length - HASHLEN;
    for (uint32_t j = 0; j < size; ++j) {
      uint32_t hash_value = Kmer2Integer(&(protein.sequence[j]));
      kmer_count[hash_value]++;
      kmer_db_exist.insert(hash_value);
    }
  }

  /* allocate memory for each kmer */
  for (std::tr1::unordered_map<uint32_t, uint32_t>::const_iterator it =
      kmer_count.begin(); it != kmer_count.end(); ++it) {
    kmer_dblocations[it->first].resize(kmer_count[it->first]);
  }

  /* Hash to Bucket */
  kmer_count.clear();
  for (uint32_t i = 0; i < proteindb.num_of_proteins; ++i) {
    const CProtein& protein = proteindb.proteins[i];
    if (protein.length < HASHLEN)
      continue;
    uint32_t size = protein.length - HASHLEN;
    for (uint32_t j = 0; j < size; ++j) {
      uint32_t hash_value = Kmer2Integer(&(protein.sequence[j]));
      kmer_dblocations[hash_value][kmer_count[hash_value]] = DBLocation(i, j);
      kmer_count[hash_value]++;
    }
  }
}

void BuildKmerNeighbors(KMERNEIGHBORS& kmer_neighbors,
                        const KMER_DB_EXIST& kmer_db_exist) {
  INFO("BUILD KMER NEIGHBORS...");
  uint32_t num_of_kmer = static_cast<uint32_t>(pow(ALPHABETSIZE, HASHLEN));
  uint32_t max_num_of_neighbors;
  Option::GetOption("-m", max_num_of_neighbors, 100);
  KNearestNeighbor k_nearest_neighbor;
  for (uint32_t i = 0; i < num_of_kmer; ++i) {
    k_nearest_neighbor.UpdateWeight(i);
    uint32_t path_id = 0, candidate = 0;
    double raw_score = 0.0;
    bool bpath_exist = true;
    while (bpath_exist && path_id < max_num_of_neighbors) {
      bpath_exist = k_nearest_neighbor.FindCandidate(candidate, path_id,
                                                     raw_score);
      if (bpath_exist && kmer_db_exist.find(candidate) != kmer_db_exist.end()) {
        kmer_neighbors[i].push_back(candidate);
      }
      path_id++;
    }
  }
}

void WriteIndex(const CProteinDB& proteindb,
                const KMERDBLOCATIONS& kmer_dblocations,
                const KMERNEIGHBORS& kmer_neighbors) {
  string index_output_file;
  Option::GetOption("-o", index_output_file);
  FILE* fout = fopen(index_output_file.c_str(), "wb");
  FILE_OPEN_CHECK(fout);
  INFO("WRITE INDEX", index_output_file.c_str());

  // proteindb
  fwrite(&(proteindb.db_length), sizeof(uint64_t), 1, fout);
  fwrite(&(proteindb.num_of_proteins), sizeof(uint32_t), 1, fout);
  fwrite(&(proteindb.max_protein_length), sizeof(uint32_t), 1, fout);
  fwrite(&(proteindb.min_protein_length), sizeof(uint32_t), 1, fout);
  for (uint32_t i = 0; i < proteindb.num_of_proteins; ++i) {
    const CProtein& protein = proteindb.proteins[i];
    uint32_t protein_name_len = protein.name.size();
    if (protein_name_len > 255) {
      protein_name_len = 255;
    }
    fwrite(&protein_name_len, sizeof(uint32_t), 1, fout);
    fwrite(protein.name.c_str(), sizeof(char), protein_name_len, fout);
    fwrite(&(protein.length), sizeof(uint32_t), 1, fout);
    fwrite(&(protein.sequence[0]), sizeof(char), protein.length, fout);
  }

  // kmer_dblocations
  uint32_t num_of_keys = kmer_dblocations.size();
  fwrite(&(num_of_keys), sizeof(uint32_t), 1, fout);
  for (KMERDBLOCATIONS::const_iterator it = kmer_dblocations.begin();
      it != kmer_dblocations.end(); ++it) {
    uint32_t hash_key = it->first;
    fwrite(&(hash_key), sizeof(uint32_t), 1, fout);
    uint32_t num_of_values = it->second.size();
    fwrite(&(num_of_values), sizeof(uint32_t), 1, fout);
    fwrite(&(it->second[0]), sizeof(DBLocation), num_of_values, fout);
  }

  // kmer_neighbors
  num_of_keys = kmer_neighbors.size();
  fwrite(&(num_of_keys), sizeof(uint32_t), 1, fout);
  for (KMERNEIGHBORS::const_iterator it = kmer_neighbors.begin();
      it != kmer_neighbors.end(); ++it) {
    uint32_t hash_key = it->first;
    fwrite(&(hash_key), sizeof(uint32_t), 1, fout);
    uint32_t num_of_values = it->second.size();
    fwrite(&(num_of_values), sizeof(uint32_t), 1, fout);
    fwrite(&(it->second[0]), sizeof(uint32_t), num_of_values, fout);
  }

  fclose(fout);
}

void ReadIndex(CProteinDB& proteindb, KMERDBLOCATIONS& kmer_dblocations,
               KMERNEIGHBORS& kmer_neighbors) {
  string index_file;
  Option::GetOption("-i", index_file);
  FILE* fin = fopen(index_file.c_str(), "rb");
  FILE_OPEN_CHECK(fin);
  INFO("READ INDEX", index_file.c_str());

  // proteindb
  FREAD_CHECK(fread(&(proteindb.db_length), sizeof(uint64_t), 1, fin), 1);
  FREAD_CHECK(fread(&(proteindb.num_of_proteins), sizeof(uint32_t), 1, fin), 1);
  FREAD_CHECK(fread(&(proteindb.max_protein_length), sizeof(uint32_t), 1, fin),
              1);
  FREAD_CHECK(fread(&(proteindb.min_protein_length), sizeof(uint32_t), 1, fin),
              1);
  proteindb.proteins.resize(proteindb.num_of_proteins);
  char protein_name[256];
  uint32_t protein_name_len;
  for (uint32_t i = 0; i < proteindb.num_of_proteins; ++i) {
    CProtein& protein = proteindb.proteins[i];
    FREAD_CHECK(fread(&protein_name_len, sizeof(uint32_t), 1, fin), 1);
    FREAD_CHECK(fread(protein_name, sizeof(char), protein_name_len, fin),
                protein_name_len);
    protein_name[protein_name_len] = 0;
    protein.name = protein_name;
    FREAD_CHECK(fread(&(protein.length), sizeof(uint32_t), 1, fin), 1);
    protein.sequence.resize(protein.length);
    FREAD_CHECK(
        fread(&(protein.sequence[0]), sizeof(char), protein.length, fin),
        protein.length);
  }

  // kmer_dblocations
  uint32_t num_of_keys = 0, num_of_values = 0, hash_key = 0;
  FREAD_CHECK(fread(&num_of_keys, sizeof(uint32_t), 1, fin), 1);
  for (uint32_t i = 0; i < num_of_keys; ++i) {
    FREAD_CHECK(fread(&hash_key, sizeof(uint32_t), 1, fin), 1);
    FREAD_CHECK(fread(&num_of_values, sizeof(uint32_t), 1, fin), 1);
    vector<DBLocation> hash_values(num_of_values);
    FREAD_CHECK(
        fread(&(hash_values[0]), sizeof(DBLocation), num_of_values, fin),
        num_of_values);
    kmer_dblocations.insert(make_pair(hash_key, hash_values));
  }

  // kmer_neighbors
  FREAD_CHECK(fread(&num_of_keys, sizeof(uint32_t), 1, fin), 1);
  for (uint32_t i = 0; i < num_of_keys; ++i) {
    FREAD_CHECK(fread(&hash_key, sizeof(uint32_t), 1, fin), 1);
    FREAD_CHECK(fread(&num_of_values, sizeof(uint32_t), 1, fin), 1);
    vector<uint32_t> hash_values(num_of_values);
    FREAD_CHECK(
        fread(&(hash_values[0]), sizeof(uint32_t), num_of_values, fin),
        num_of_values);
    kmer_neighbors.insert(make_pair(hash_key, hash_values));
  }

  fclose(fin);
}
