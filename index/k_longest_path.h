#ifndef KLONGESTPATH_H_
#define KLONGESTPATH_H_

#include "./../util/bio_util.h"
#include <queue>
#include <algorithm>

namespace longest_path {

#define NODE_NUM 200
#define MAX_K  5001

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

class KLongestPath {
 private:
  priority_queue<QueNodeDP> cand[NODE_NUM];
  int *dp[NODE_NUM], *pre[NODE_NUM];
  struct AdjLinkDP adj[100010], invAdj[100010];
  int adjN;
  int invAdjN;

  int nNumNode;
  vector<int> vNodeLabel;
  vector<vector<usint32_t> > vLinkList;
  vector<usint32_t> vLinkNodeSize;
  vector<vector<int> > vEdgeWeight;

 public:
  KLongestPath() {
    for (usint32_t i = 0; i < NODE_NUM; ++i) {
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
    for (usint32_t i = 0; i < vLinkList.size(); i++) {
      vLinkList[i].resize(70);
      vEdgeWeight[i].resize(70);
    }
    BuildDAG();
  }
  ~KLongestPath() {
    INFO("Release K longest path memory...");
    for (usint32_t i = 0; i < NODE_NUM; ++i) {
      delete[] dp[i];
      delete[] pre[i];
    }
  }

 private:
  void insert(AdjLinkDP * adj, int &adjN, int a, int b, int c);
  int Query(int u, int k);
  void ShowPath(int n, int k, vector<usint32_t> & pathTmp);
  int GetWeight(const char & queryAA, const char & AA);
  void BuildDAG();

 public:
  void UpdateWeight(const char * querySeq);
  bool FindCandidate(char * strCandPep, const int & i, int & rawScore);
};

}  // namespace longest_path
#endif /* KLONGESTPATH_H_ */
