#pragma once
#ifndef BUILDINDEX_H_
#define BUILDINDEX_H_

#include "sdk.h"
#include "option.h"
#include "KlongestPathDP.h"
#include <map>
#include <set>
#include <vector>
#include <string>

typedef struct {
	vector<string> vProName;
	vector<usint32_t> vProSize;
	char** vProSeq;
	usint32_t nNumPro;
	usint64_t sizeAA;
} CProteinDB;

typedef struct {
	usint64_t * counter;
	usint32_t * index;
	usint32_t nSizeCounter;
	usint64_t nSizeIndex;
} CKmerNearest;

typedef struct {
	usint64_t * counter;
	usint32_t * index_proID;
	usint32_t * index_proPos;
	usint32_t nSizeCounter;
	usint64_t nSizeIndex;
} CKmerPosition;

class Index {
	public:
		Index(const Option& opt);
		~Index();

		void BuildIndex();
		void ReadIndex();

	private:
		void AnalyzeProteinDataBase(const char * strRef, const usint64_t & refLen);
		void BuildNearestKmer();
		void BuildKmerPosition();

		void CountNearestKmer();
		void CountKmerPosition();

		void WriteIndex();

		void GetThresholdScore();

		void TestIndex();

	public:
		//unordered_map<usint32_t, vector<usint32_t> > kmer_nearest_;
		//unordered_map<usint32_t, vector<pair<usint32_t, usint32_t> > > kmer_position_;
		CKmerNearest kmer_nearest_;
		CKmerPosition kmer_position_;
		CProteinDB proDB_;
	private:
		Option opt_;
		CKlongestPathDP m_kLongestPath_;

		int score_threshold_;
};

#endif /* BUILDINDEX_H_ */
