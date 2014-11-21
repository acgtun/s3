#pragma once
#ifndef DBSEARCH_H_
#define DBSEARCH_H_

#include "sdk.h"
#include "option.h"
#include "buildindex.h"

typedef struct {
	char ** readsArray;
	vector<usint32_t> vReadSize;
	usint32_t nReadsNum;
} CReads;

class DBSearch {
public:
	DBSearch(const Option& opt, const Index * index);
	~DBSearch();

	void Search();

private:
	void ReadReads();

private:
	Option opt_;
	CReads reads_;
	usint32_t seed_start_pos_[MAX_QUERY_LEN / 6 + 1];
	const Index * index_;
};

#endif /* DBSEARCH_H_ */
