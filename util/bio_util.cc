#include "bio_util.h"

const usint32_t base_preCal[20][6] = {
    { 0 * basep[0], 0 * basep[1], 0 * basep[2], 0 * basep[3], 0 * basep[4], 0 * basep[5] },
    { 1 * basep[0], 1 * basep[1], 1 * basep[2], 1 * basep[3], 1 * basep[4], 1 * basep[5] },
    { 2 * basep[0], 2 * basep[1], 2 * basep[2], 2 * basep[3], 2 * basep[4], 2 * basep[5] },
    { 3 * basep[0], 3 * basep[1], 3 * basep[2], 3 * basep[3], 3 * basep[4], 3 * basep[5] },
    { 4 * basep[0], 4 * basep[1], 4 * basep[2], 4 * basep[3], 4 * basep[4], 4 * basep[5] },
    { 5 * basep[0], 5 * basep[1], 5 * basep[2], 5 * basep[3], 5 * basep[4], 5 * basep[5] },
    { 6 * basep[0], 6 * basep[1], 6 * basep[2], 6 * basep[3], 6 * basep[4], 6 * basep[5] },
    { 7 * basep[0], 7 * basep[1], 7 * basep[2], 7 * basep[3], 7 * basep[4], 7 * basep[5] },
    { 8 * basep[0], 8 * basep[1], 8 * basep[2], 8 * basep[3], 8 * basep[4], 8 * basep[5] },
    { 9 * basep[0], 9 * basep[1], 9 * basep[2], 9 * basep[3], 9 * basep[4], 9 * basep[5] },
    {10 * basep[0],10 * basep[1],10 * basep[2],10 * basep[3],10 * basep[4],10 * basep[5] },
    {11 * basep[0],11 * basep[1],11 * basep[2],11 * basep[3],11 * basep[4],11 * basep[5] },
    {12 * basep[0],12 * basep[1],12 * basep[2],12 * basep[3],12 * basep[4],12 * basep[5] },
    {13 * basep[0],13 * basep[1],13 * basep[2],13 * basep[3],13 * basep[4],13 * basep[5] },
    {14 * basep[0],14 * basep[1],14 * basep[2],14 * basep[3],14 * basep[4],14 * basep[5] },
    {15 * basep[0],15 * basep[1],15 * basep[2],15 * basep[3],15 * basep[4],15 * basep[5] },
    {16 * basep[0],16 * basep[1],16 * basep[2],16 * basep[3],16 * basep[4],16 * basep[5] },
    {17 * basep[0],17 * basep[1],17 * basep[2],17 * basep[3],17 * basep[4],17 * basep[5] },
    {18 * basep[0],18 * basep[1],18 * basep[2],18 * basep[3],18 * basep[4],18 * basep[5] },
    {19 * basep[0],19 * basep[1],19 * basep[2],19 * basep[3],19 * basep[4],19 * basep[5] }
};

void DeCodeAA(const usint32_t & hashValue, char * strSeed) {
  usint32_t value = hashValue;
  for (usint32_t i = 0; i < HASHAALEN; i++) {
    strSeed[i] = AA20[value % 20];
    value /= 20;
  }
  strSeed[HASHAALEN] = 0;
}

usint32_t GetHashValue(const char * strSeed) {
  usint32_t hashValue = 0;
  for (usint32_t i = 0; i < HASHAALEN; i++) {
    hashValue += base_preCal[base[strSeed[i] - 'A']][i];
  }
  return hashValue;
}

usint32_t GetPower(const usint32_t & base, const usint32_t & p) {
  usint32_t ret = 1;
  for (usint32_t i = 0; i < p; i++) {
    ret *= base;
  }
  return ret;
}

void DisplayResults(const char* query_name, const char* database_file,
                    const vector<M8Results>& aligned_results, const int& outfmt,
                    ofstream& fout) {
  if (outfmt == 7) {
    fout << "# RMP 1.0.0 Oct, 2014" << endl;
    fout << "# Query: " << query_name << endl;
    fout << "# Database: " << database_file << endl;
    fout
        << "# Fields: query id, subject id, % identity, alignment length, mismatches, "
            "gap opens, q. start, q. end, s. start, s. end, evalue, bit score"
        << endl;
    fout << "# " << aligned_results.size() << " hits found" << endl;
  }
  for (usint32_t i = 0; i < aligned_results.size(); i++) {
    fout << query_name << "\t" << aligned_results[i].protein_name << "\t"
         << aligned_results[i].identity << "\t"
         << aligned_results[i].aligned_len << "\t"
         << aligned_results[i].mismatch << "\t" << aligned_results[i].gap_open
         << "\t" << aligned_results[i].qs << "\t" << aligned_results[i].qe
         << "\t" << aligned_results[i].ps << "\t" << aligned_results[i].pe
         << "\t" << aligned_results[i].evalue << "\t"
         << aligned_results[i].bit_score << endl;
  }
}
