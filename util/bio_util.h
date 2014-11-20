#ifndef BIO_UTIL_H_
#define BIO_UTIL_H_

#include "sdk.h"

#include <vector>
#include <fstream>
#include <iostream>

using namespace std;

/*
 A0
 R1
 N2
 D3
 C4
 E5
 Q6
 G7
 H8
 I9
 L10
 K11
 M12
 F13
 P14
 S15
 T16
 W17
 Y18
 V19
 */
/* B J O U X Z*/

#define HASHAALEN 6 
const string AA20 = "ARNDCEQGHILKMFPSTWYV";
const string NODELABEL = " ARNDCEQGHILKMFPSTWYV";
                   /*A   B  C  D  E   F  G  H  I   J   K   L   M  N   O   P  Q  R   S   T   U   V   W   X   Y   Z */
const int base[] = { 0, -1, 4, 3, 6, 13, 7, 8, 9, -1, 11, 10, 12, 2, -1, 14, 5, 1, 15, 16, -1, 19, 17, -1, 18, -1 };
const usint32_t basep[] = { 1, 20, 400, 8000, 160000, 3200000, 64000000, 1280000000 };

const int BLOSUM62[][20] = {
//A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V
{ 4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0 },  //A
{ -1, 5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3 },  //R
{ -2, 0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3 },  //N
{ -2,-2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3 },  //D
{  0,-3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1 },  //C
{ -1, 1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2 },  //Q
{ -1, 0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2 },  //E
{ 0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3 },  //G
{ -2, 0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3 },  //H
{ -1,-3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3 },  //I
{ -1,-2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1 },  //L
{ -1, 2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2 },  //K
{ -1,-1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1 },  //M
{ -2,-3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1 },  //F
{ -1,-2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2 },  //P
{ 1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2 },  //S
{ 0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0 },  //T
{ -3,-3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3 },  //W
{ -2,-2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1 },  //Y
{ 0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3,  -1, 4 } };  //V
//  http://www.ncbi.nlm.nih.gov/BLAST/blastcgihelp.shtml#get_subsequence
//  http://www.ncbi.nlm.nih.gov/Class/FieldGuide/BLOSUM62.txt

/* Given a 6-mer amino acids sequence, translate it to a integer number.
 * Use 20-based number, A is 0, R is 1, N is 3 and so on.
 * */
usint32_t GetHashValue(const char * strSeed);

/* Given a integer, translate it to a 6-mer amino acids, it's also 20-based.*/
void DeCodeAA(const usint32_t & hashValue, char * strSeed);

/* calculate base^p */
usint32_t GetPower(const usint32_t & base, const usint32_t & p);

/* M8Results is a data structure to store the results of a protein alingment, same as BLAST -m 8*/
struct M8Results {
  M8Results(const string& pro_name, const double& idty, const int& ali_len,
            const int& mis, const int& gap, const usint32_t& q_start,
            const usint32_t& q_end, const usint32_t& p_start,
            const usint32_t& p_end, const double& e_value,
            const double& bitScore)
      : protein_name(pro_name),
        identity(idty),
        aligned_len(ali_len),
        mismatch(mis),
        gap_open(gap),
        qs(q_start),
        qe(q_end),
        ps(p_start),
        pe(p_end),
        evalue(e_value),
        bit_score(bitScore) {
  }
  M8Results() {
    protein_name = "";

    identity = 0.0;
    aligned_len = 0;
    mismatch = 0;
    gap_open = 0;

    qs = 0;
    qe = 0;
    ps = 0;
    pe = 0;

    evalue = 0.0;
    bit_score = 0.0;
  }

  string protein_name;

  double identity;
  int aligned_len;
  int mismatch;
  int gap_open;

  int qs;
  int qe;
  int ps;
  int pe;

  double evalue;
  double bit_score;

  static bool SORT_CMP_EValue(const M8Results& a, const M8Results& b) {
    return a.evalue < b.evalue;
  }
};

/* Output the alingment results for one query to fout */
void DisplayResults(const char* query_name, const char* database_file,
                    const vector<M8Results>& aligned_results, const int& outfmt,
                    ofstream& fout);

#endif /* BIO_UTIL_H_ */
