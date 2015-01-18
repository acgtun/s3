#ifndef BIO_UTIL_H_
#define BIO_UTIL_H_

#include "sdk.h"

#include <string>
#include <iostream>

using std::cout;
using std::cerr;
using std::endl;
using std::string;


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

#define HASHLEN 6
#define ALPHABETSIZE 8

const string AA20 = "ARNDCEQGHILKMFPSTWYV";
const string NODELABEL = " ARNDCEQGHILKMFPSTWYV";

const int AAINDEX[] =
{ 0, -1, 4, 3, 6, 13, 7, 8, 9, -1, 11, 10, 12, 2, -1, 14, 5, 1, 15, 16, -1, 19, 17, -1, 18, -1 };
/*A   B  C  D  E   F  G  H  I   J   K   L   M  N   O   P  Q  R   S   T   U   V   W   X   Y   Z */

// [A S T] [R K E D Q] [N H] [C] [G] [I V L M] [F Y W] [P]
//    0         1        2    3   4      5        6     7
const int REDUCEDAAINDEX[] =
{ 0, -1, 3, 1, 1, 6, 4, 2, 5, -1, 1, 5, 5, 2, -1, 7, 1, 1, 0, 0, -1, 5, 6, -1, 6, -1 };
/*A   B  C  D  E  F  G  H  I   J  K  L  M  N   O  P  Q  R  S  T   U  V  W   X  Y   Z */

const uint32_t BASEP[] = { 1, 8, 64, 512, 4096, 32768, 262144, 2097152 };

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

const double REDUCEDBLOSUM62[ALPHABETSIZE][ALPHABETSIZE] = {
{ 1.88889,  -0.8,     -1,       -0.666667, -0.666667, -1.08333, -2.22222, -1       },
{-0.8,       1.52,    -0.1,     -3.2,      -1.8,      -2.35,    -2.66667, -1.2     },
{-1,        -0.1,      4,       -3,        -1,        -2.75,    -1.66667, -2       },
{-0.666667, -3.2,     -3,        9,        -3,        -1,       -2,       -3       },
{-0.666667, -1.8,     -1,       -3,         6,        -3.5,     -2.66667, -2       },
{-1.08333,  -2.35,    -2.75,    -1,        -3.5,       2.3125,  -1.16667, -2.5     },
{-2.22222,  -2.66667, -1.66667, -2,        -2.66667,  -1.16667,  4,       -3.66667 },
{-1,        -1.2,     -2,       -3,        -2,        -2.5,     -3.66667,  7       } };

inline uint32_t Kmer2Integer(const char* kmer) {
  uint32_t hash_value = 0;
  for (uint32_t i = 0; i < HASHLEN; ++i) {
    hash_value += REDUCEDAAINDEX[kmer[i] - 'A'] * BASEP[i];
  }
  return hash_value;
}

/* transfer the integer to 8-based number */
inline string Integer2KmerDigit(const uint32_t& hash_value) {
  string kmer;
  uint32_t n = hash_value, j = 0;
  while (n) {
    kmer += 48 + n % ALPHABETSIZE;
    j++;
    n /= ALPHABETSIZE;
  }
  while (j < HASHLEN) {
    kmer += 48;
    j++;
  }
  return kmer;
}

#endif /* BIO_UTIL_H_ */
