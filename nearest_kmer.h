#ifndef NEAREST_KMER_H_
#define NEAREST_KMER_H_

#include "bio_util.h"

#include <vector>
#include <queue>
#include <algorithm>

using std::vector;
using std::priority_queue;

namespace nearest_kmer {

#define NODE_NUM 200
#define MAX_K  5001

class AdjLinkDP {
 public:
  int v, nxt;
  double c;
};

class QueNodeDP {
 public:
  int pre;
  double c, tc;
  bool operator<(const QueNodeDP & b) const {
    return c < b.c;
  }
};

class KNearestNeighbor {
 private:
  priority_queue<QueNodeDP> cand[NODE_NUM];
  double *dp[NODE_NUM];
  int *pre[NODE_NUM];
  int dp0[NODE_NUM];
  struct AdjLinkDP adj[100010], invAdj[100010];
  int adjN;
  int invAdjN;

  int nNumNode;
  vector<int> vNodeLabel;
  vector<vector<uint32_t> > vLinkList;
  vector<uint32_t> vLinkNodeSize;
  vector<vector<double> > vEdgeWeight;

 public:
  KNearestNeighbor() {
    for (uint32_t i = 0; i < NODE_NUM; ++i) {
      dp[i] = new double[MAX_K];
      pre[i] = new int[MAX_K];
    }
    adjN = 0;
    invAdjN = 0;
    nNumNode = 0;

    vNodeLabel.resize(NODE_NUM);
    vLinkList.resize(NODE_NUM);
    vEdgeWeight.resize(NODE_NUM);
    vLinkNodeSize.resize(NODE_NUM);
    for (uint32_t i = 0; i < vLinkList.size(); i++) {
      vLinkList[i].resize(70);
      vEdgeWeight[i].resize(70);
    }
    BuildDAG();
  }
  ~KNearestNeighbor() {
    for (uint32_t i = 0; i < NODE_NUM; ++i) {
      delete[] dp[i];
      delete[] pre[i];
    }
  }

 private:
  void insert(AdjLinkDP* adj, int& adjN, const int& a, const int& b,
              const double& c);
  double Query(int u, int k);
  void ShowPath(int n, int k, vector<uint32_t> & pathTmp);
  void BuildDAG();

 public:
  void UpdateWeight(const uint32_t& kmer);
  bool FindCandidate(uint32_t& candidate, const int& i, double& rawScore);
};

}  // namespace nearest_kmer
#endif /* NEAREST_KMER_H_ */
