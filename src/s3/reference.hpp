/*
 *    This is the header file for building index.
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

#ifndef REFERENCE_H_
#define REFERENCE_H_

#include "util.hpp"

struct CProtein {
  string name;
  uint32_t length;
  string sequence;
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
void BuildKmerLocation(const CProteinDB& proteindb, HashTable& hash_table);

void WriteIndex(const CProteinDB& proteindb, const HashTable& hash_table,
                const string& output_index_file);
void ReadIndex(CProteinDB& proteindb, HashTable& hash_table);

#endif /* REFERENCE_H_ */
