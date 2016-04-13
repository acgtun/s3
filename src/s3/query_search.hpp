/*
 *    This is the header file of searching proteins.
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

#ifndef QUERY_SEARCH_H_
#define QUERYS_EARCH_H_

#include "reference.hpp"
#include "util.hpp"

#include <vector>
#include <fstream>

using std::vector;
using std::ofstream;

void GetTopProteinIDS(const HashTable& hash_table,
                      const char* query_seq,
                      const uint32_t& num_of_top_proteins,
                      vector<uint32_t>& protein_candidates_id,
                      uint32_t& num_of_protein_candidates);

#endif /* QUERY_SEARCH_H_ */
