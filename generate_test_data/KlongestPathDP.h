/*
 * DPKlongestPath.h
 *
 *  Created on: 2011-8-15
 *      Author: chenhaifeng
 */
#ifndef DPKLONGESTPATH_H_
#define DPKLONGESTPATH_H_

#include "sdk.h"
#include <string>
#include <stdio.h>
#include <queue>
#include <time.h>
#include <fstream>
#include <iostream>
using namespace std;

#define NODE_NUM 200
#define MAX_K  5001

typedef struct {
	char chr[3];
} NODElABEL;

class AdjLinkDP {
public:
	int v, c, nxt;
};

class QueNodeDP {
public:
	int c, tc, pre;
	bool operator<(const QueNodeDP & b) const {
		return c < b.c;
	}
};

class CKlongestPathDP {
private:
	priority_queue<QueNodeDP> cand[NODE_NUM];
	int *dp[NODE_NUM], *pre[NODE_NUM];
	struct AdjLinkDP adj[100010], invAdj[100010];
	int adjN;
	int invAdjN;

	int nNumNode;
	vector<int> vNodeLabel;
	vector<vector<size_t> > vLinkList;
	vector<size_t> vLinkNodeSize;
	vector<vector<int> > vEdgeWeight;

public:
	CKlongestPathDP() {
		for (size_t i = 0; i < NODE_NUM; ++i) {
			dp[i] = new int[MAX_K];
			pre[i] = new int[MAX_K];
		}
		adjN = 0;
		invAdjN = 0;
		nNumNode = 0;

		vNodeLabel.resize(NODE_NUM);
		vLinkList.resize(NODE_NUM);
		vEdgeWeight.resize(NODE_NUM);
		vLinkNodeSize.resize(NODE_NUM);
		for (size_t i = 0; i < vLinkList.size(); i++) {
			vLinkList[i].resize(70);
			vEdgeWeight[i].resize(70);
		}
		BuildDAG();
	}
	~CKlongestPathDP() {
		for (size_t i = 0; i < NODE_NUM; ++i) {
			delete[] dp[i];
			delete[] pre[i];
		}
	}

private:
	void insert(AdjLinkDP * adj, int &adjN, int a, int b, int c);
	int Query(int u, int k);
	void ShowPath(int n, int k, vector<size_t> & pathTmp);
	int GetWeight(const char & queryAA, const char & AA);
	void BuildDAG();

public:

	void UpdateWeight(const char * querySeq);
	bool FindCandidiate(char * strCandPep, const int & i, int & rawScore);
	//void DPKlongestPath(const char * querySeq, const int & thresholdScore, vector<string> & vCandPep);
};
#endif /* DPKLONGESTPATH_H_ */
