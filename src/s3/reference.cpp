/*
 *    This file contains functions for building index.
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

#include "util.hpp"
#include "reference.hpp"

void ReadFASTAFile(const string& fasta_file, vector<string>& names,
                   vector<string>& seqs) {
  fprintf(stderr, "[READ FASTA FILE %s]\n", fasta_file.c_str());
  ifstream fin(fasta_file.c_str());
  if (!fin.good()) {
    throw SMITHLABException("cannot open the file " + fasta_file);
  }
  string sequence, line;
  while (getline(fin, line)) {
    if (line[0] == '>') {
      if (sequence.size() != 0) {
        seqs.push_back(sequence);
        sequence.clear();
      }
      uint32_t space_pos = line.find_first_of(' ');
      if (space_pos == string::npos) {
        names.push_back(line.substr(1));
      } else {
        names.push_back(line.substr(1, space_pos - 1));
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
    seqs.push_back(sequence);
  }
  fin.close();
}

void BuildProteinDB(const string& database_file, CProteinDB& proteindb) {
  vector<string> protein_names, protein_seqs;
  ReadFASTAFile(database_file, protein_names, protein_seqs);

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
  fprintf(stderr, "THE NUMBER OF PROTEINS IN DATABASE IS %u\n",
          proteindb.num_of_proteins);
  fprintf(stderr, "THE NUMBER OF AMINO ACIDES IN DATABASE IS %lu\n",
          proteindb.db_length);
}

void BuildKmerLocation(const CProteinDB& proteindb, HashTable& hash_table) {
  fprintf(stderr, "BUILD KMER LOCATIONS...");
  for (uint32_t i = 0; i < proteindb.num_of_proteins; ++i) {
    const CProtein& protein = proteindb.proteins[i];
    if (protein.length < HASHLEN)
      continue;
    uint32_t size = protein.length - HASHLEN;
    for (uint32_t j = 0; j < size; ++j) {
      uint32_t hash_value = Kmer2Integer(&(protein.sequence[j]));
      hash_table[hash_value].push_back(make_pair(i, j));
    }
  }
}

void WriteIndex(const CProteinDB& proteindb, const HashTable& hash_table,
                const string& outfile) {
  fprintf(stderr, "WRITE INDEX TO %s\n", outfile.c_str());
  FILE* fout = fopen(outfile.c_str(), "wb");
  if (!fout) {
    throw SMITHLABException("cannot open input file " + outfile);
  }

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

  // hash_table
  uint32_t num_of_keys = hash_table.size();
  fwrite(&(num_of_keys), sizeof(uint32_t), 1, fout);
  for (HashTable::const_iterator it = hash_table.begin();
      it != hash_table.end(); ++it) {
    uint32_t hash_key = it->first;
    fwrite(&(hash_key), sizeof(uint32_t), 1, fout);
    uint32_t num_of_values = it->second.size();
    fwrite(&(num_of_values), sizeof(uint32_t), 1, fout);
    fwrite(&(it->second[0]), sizeof(pair<uint32_t, uint32_t> ), num_of_values,
           fout);
  }

  fclose(fout);
}

void ReadIndex(CProteinDB& proteindb, HashTable& hash_table,
               const string& index_file) {
  fprintf(stderr, "READ INDEX FROM %s\n", index_file.c_str());
  FILE* fin = fopen(index_file.c_str(), "rb");
  if (!fin) {
    throw SMITHLABException("cannot open input file " + index_file);
  }

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
  fprintf(stderr, "THE NUMBER OF PROTEINS IN DATABASE IS %u\n",
          proteindb.num_of_proteins);
  fprintf(stderr, "THE NUMBER OF AMINO ACIDES IN DATABASE IS %lu\n",
          proteindb.db_length);

  // hash_table
  uint32_t num_of_keys = 0, num_of_values = 0, hash_key = 0;
  FREAD_CHECK(fread(&num_of_keys, sizeof(uint32_t), 1, fin), 1);
  for (uint32_t i = 0; i < num_of_keys; ++i) {
    FREAD_CHECK(fread(&hash_key, sizeof(uint32_t), 1, fin), 1);
    FREAD_CHECK(fread(&num_of_values, sizeof(uint32_t), 1, fin), 1);
    vector < pair<uint32_t, uint32_t> > hash_values(num_of_values);
    FREAD_CHECK(
        fread(&(hash_values[0]), sizeof(pair<uint32_t, uint32_t> ),
              num_of_values, fin),
        num_of_values);
    hash_table.insert(make_pair(hash_key, hash_values));
  }

  fclose(fin);
}
