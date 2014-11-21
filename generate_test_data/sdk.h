#pragma once
#ifndef SDK_H_
#define SDK_H_

#include <time.h>
#include <math.h>
#include <ctype.h>
#include <stdio.h>
#include <errno.h>
#include <limits.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>

#include <map>
#include <set>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <algorithm>

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

typedef uint64_t usint64_t;
typedef uint32_t usint32_t;
#define MAX_LINE_LEN 1024
#define HASHAALEN 6
#define MAX_QUERY_LEN 10000


const double GB = 1024.00 * 1024.00 * 1024.00;
const string AA20 = "ARNDCEQGHILKMFPSTWYV";
const string NODELABEL = " ARNDCEQGHILKMFPSTWYV";

const int base[] = { 0, -1, 4, 3, 5, 13, 7, 8, 9, -1, 11, 10, 12, 2, -1, 14, 6, 1, 15, 16, -1, 19, 17, -1, 19, -1 };
const usint32_t basep[] = { 1, 20, 400, 8000, 160000, 3200000, 64000000, 1280000000 };
const usint32_t base_preCal[20][6] = {
{ 0 * basep[0],  0 * basep[1],  0 * basep[2],  0 *  basep[3],  0 * basep[4],  0 * basep[5]},
{ 1 * basep[0],  1 * basep[1],  1 * basep[2],  1 *  basep[3],  1 * basep[4],  1 * basep[5]},
{ 2 * basep[0],  2 * basep[1],  2 * basep[2],  2 *  basep[3],  2 * basep[4],  2 * basep[5]},
{ 3 * basep[0],  3 * basep[1],  3 * basep[2],  3 *  basep[3],  3 * basep[4],  3 * basep[5]},
{ 4 * basep[0],  4 * basep[1],  4 * basep[2],  4 *  basep[3],  4 * basep[4],  4 * basep[5]},
{ 5 * basep[0],  5 * basep[1],  5 * basep[2],  5 *  basep[3],  5 * basep[4],  5 * basep[5]},
{ 6 * basep[0],  6 * basep[1],  6 * basep[2],  6 *  basep[3],  6 * basep[4],  6 * basep[5]},
{ 7 * basep[0],  7 * basep[1],  7 * basep[2],  7 *  basep[3],  7 * basep[4],  7 * basep[5]},
{ 8 * basep[0],  8 * basep[1],  8 * basep[2],  8 *  basep[3],  8 * basep[4],  8 * basep[5]},
{ 9 * basep[0],  9 * basep[1],  9 * basep[2],  9 *  basep[3],  9 * basep[4],  9 * basep[5]},
{10 * basep[0], 10 * basep[1], 10 * basep[2], 10 *  basep[3], 10 * basep[4], 10 * basep[5]},
{11 * basep[0], 11 * basep[1], 11 * basep[2], 11 *  basep[3], 11 * basep[4], 11 * basep[5]},
{12 * basep[0], 12 * basep[1], 12 * basep[2], 12 *  basep[3], 12 * basep[4], 12 * basep[5]},
{13 * basep[0], 13 * basep[1], 13 * basep[2], 13 *  basep[3], 13 * basep[4], 13 * basep[5]},
{14 * basep[0], 14 * basep[1], 14 * basep[2], 14 *  basep[3], 14 * basep[4], 14 * basep[5]},
{15 * basep[0], 15 * basep[1], 15 * basep[2], 15 *  basep[3], 15 * basep[4], 15 * basep[5]},
{16 * basep[0], 16 * basep[1], 16 * basep[2], 16 *  basep[3], 16 * basep[4], 16 * basep[5]},
{17 * basep[0], 17 * basep[1], 17 * basep[2], 17 *  basep[3], 17 * basep[4], 17 * basep[5]},
{18 * basep[0], 18 * basep[1], 18 * basep[2], 18 *  basep[3], 18 * basep[4], 18 * basep[5]},
{19 * basep[0], 19 * basep[1], 19 * basep[2], 19 *  basep[3], 19 * basep[4], 19 * basep[5]}};

const int BLOSUM62[][20] = {
//A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V
{ 4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0}, //A
{-1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3}, //R
{-2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3}, //N
{-2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3}, //D
{ 0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1}, //C
{-1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2}, //Q
{-1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2}, //E
{ 0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3}, //G
{-2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3}, //H
{-1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3}, //I
{-1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1}, //L
{-1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2}, //K
{-1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1}, //M
{-2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1}, //F
{-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2}, //P
{ 1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2}, //S
{ 0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0}, //T
{-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3}, //W
{-2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1}, //Y
{ 0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4}}; //V
//  http://www.ncbi.nlm.nih.gov/BLAST/blastcgihelp.shtml#get_subsequence
//  http://www.ncbi.nlm.nih.gov/Class/FieldGuide/BLOSUM62.txt

inline void MemoryAllocateCheck(void * pointer, const char * file, int line) {
	if (pointer == NULL) {
		printf("Memory allocate error in %s at line %d\n", file, line);
		exit(EXIT_FAILURE);
	}
}

inline void FileOpenCheck(FILE * pfile, const char * file, int line) {
	if (pfile == NULL) {
		printf("File open error in %s at line %d\n", file, line);
		exit(EXIT_FAILURE);
	}
}

#define HANDLE_ERROR(err) (HandleError( err, __FILE__, __LINE__ ))
#define LOG_INFO printf("--- %s:%s:%d\n",  __FILE__, __func__, __LINE__)
#define FILE_OPEN_CHECK(pfile) (FileOpenCheck( pfile, __FILE__, __LINE__))
#define CHECK_READ_LEN(len, nReadsNum) (checkReadLen(len, nReadsNum, __FILE__, __LINE__))
#define MEMORY_ALLOCATE_CHECK(pointer)  (MemoryAllocateCheck(pointer, __FILE__, __LINE__))
#define LOG_INFO_CPP cout << "-----" << __FILE__ << " " << __func__ << " " << __LINE__ << endl

#define FREAD_CHECK(func, size) { \
	size_t s = func; \
	if(s != size) { \
		printf("read file error. --- %s:%s:%d\n", __FILE__, __func__, __LINE__); \
		exit(EXIT_FAILURE); \
	} \
}

#define ERROR_INFO(msg) { \
	printf("--ERROR INFO-- %s --- %s:%s:%d\n", msg, __FILE__, __func__, __LINE__); \
	exit(EXIT_FAILURE); \
}

#define TIME_INFO(func, msg) { \
	clock_t start_t, end_t; \
	start_t = clock(); \
	func; \
	end_t = clock(); \
	printf("--INFO-- %s takes %.3lf seconds.\n", msg, (double) ((end_t - start_t) / CLOCKS_PER_SEC )); \
}

inline void INFO(const string & msg) {
	cout << "--INFO-- " << msg << "." << endl;
}
inline void INFO(const string & msg, const string & val) {
	cout << "--INFO-- " << msg << " " << val << "." << endl;
}
inline void INFO(const string & msg, const usint64_t &val) {
	cout << "--INFO-- " << msg << " " << val << "." << endl;
}
inline void INFO(const string & msg, const usint32_t &val) {
	cout << "--INFO-- " << msg << " " << val << "." << endl;
}
inline void INFO(const string & msg, const int &val) {
	cout << "--INFO-- " << msg << " " << val << "." << endl;
}
inline void INFO(const string & msg, const double &val) {
	cout << "--INFO-- " << msg << " " << val << "." << endl;
}
inline void INFO(const string & msg, const usint32_t & val1, const string & val2) {
	cout << "--INFO-- " << msg << " " << val1 << " " << val2 << "." << endl;
}
inline void INFO(const string & msg, const usint64_t & val1, const string & val2) {
	cout << "--INFO-- " << msg << " " << val1 << " " << val2 << "." << endl;
}

inline void DeCodeAA(const usint32_t & hashValue, char * strSeed) {
	usint32_t value = hashValue;
	for (usint32_t i = 0; i < HASHAALEN; i++) {
		strSeed[i] = AA20[value % 20];
		value /= 20;
	}
	strSeed[HASHAALEN] = 0;
}

inline usint32_t GetHashValue(const char * strSeed) {
	usint32_t hashValue = 0;
	for (usint32_t i = 0; i < HASHAALEN; i++) {
		hashValue += base_preCal[base[strSeed[i] - 'A']][i];
	}
	return hashValue;
}

inline usint32_t GetPower(const usint32_t & base, const usint32_t & p) {
	usint32_t ret = 1;
	for (usint32_t i = 0; i < p; i++) {
		ret *= base;
	}
	return ret;
}

#endif /* SDK_H_ */
