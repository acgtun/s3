#include "nearest_kmer.h"
#include <limits.h>

namespace nearest_kmer {

void KNearestNeighbor::insert(AdjLinkDP* adj, int& adjN, const int& a,
                              const int& b, const double& c) {
  adj[adjN].nxt = adj[a].nxt;
  adj[adjN].c = c;
  adj[adjN].v = b;
  adj[a].nxt = adjN;
  adjN++;
}

double KNearestNeighbor::Query(int u, int k) {
  if (dp0[u] >= k)
    return dp[u][k];
  if (dp0[u] == 0) {
    if (invAdj[u].nxt == -1) {
      if (u != 0)
        return INT_MIN;
      QueNodeDP tmp;
      tmp.c = 0;
      tmp.pre = -1;
      tmp.tc = 0;
      cand[u].push(tmp);
      dp0[u]++;
      dp[u][1] = cand[u].top().c;
      pre[u][0]++;
      pre[u][1] = cand[u].top().pre;
      return dp[u][1];
    } else {
      for (int i = invAdj[u].nxt; i != -1; i = invAdj[i].nxt) {
        int v = invAdj[i].v;
        QueNodeDP tmp;
        double res = Query(v, 1);
        if (fabs(res - INT_MIN) < 1e-6)
          continue;
        tmp.c = res + invAdj[i].c;
        tmp.pre = v + 1 * 10000;
        tmp.tc = invAdj[i].c;
        cand[u].push(tmp);
      }
      if (cand[u].size() == 0)
        return -1;
      dp0[u]++;
      dp[u][1] = cand[u].top().c;
      pre[u][0]++;
      pre[u][1] = cand[u].top().pre;
      return dp[u][1];
    }
  } else {
    if (cand[u].empty() || invAdj[u].nxt == -1)
      return INT_MIN;
    QueNodeDP cur = cand[u].top(), tmp;
    cand[u].pop();
    double res = Query(cur.pre % 10000, cur.pre / 10000 + 1);
    if (fabs(res - INT_MIN) > 1e-6) {
      tmp.c = res + cur.tc;
      tmp.pre = cur.pre % 10000 + (cur.pre / 10000 + 1) * 10000;
      tmp.tc = cur.tc;
      cand[u].push(tmp);
    }
    if (!cand[u].empty()) {
      dp0[u]++;
      dp[u][dp0[u]] = cand[u].top().c;
      pre[u][0]++;
      pre[u][pre[u][0]] = cand[u].top().pre;
      return dp[u][dp0[u]];
    } else
      return INT_MIN;
  }
}

void KNearestNeighbor::ShowPath(int n, int k, vector<uint32_t> & pathTmp) {
  if (pre[n][k] == -1) {
    pathTmp.push_back(n);
    return;
  }

  ShowPath(pre[n][k] % 10000, pre[n][k] / 10000, pathTmp);
  pathTmp.push_back(n);
}

bool KNearestNeighbor::FindCandidate(uint32_t& candidate, const int& i,
                                     double& rawScore) {
  rawScore = Query(nNumNode - 1, i + 1);
  if (fabs(rawScore - INT_MIN) > 1e-6) {
    vector <uint32_t> pathTmp;
    ShowPath(nNumNode - 1, i + 1, pathTmp);
    candidate = 0;
    //cout << endl;
    //cout << "cand=";
    for (uint32_t i = 1; i < pathTmp.size() - 1; i++) {
      //cout << vNodeLabel[pathTmp[i]] - 1;
      candidate += (vNodeLabel[pathTmp[i]] - 1) * BASEP[i];
    }
   // cout << endl;
    return true;
  }
  return false;
}

void KNearestNeighbor::UpdateWeight(const uint32_t& kmer) {
  /* set edges' weight */
  //cout << endl;
  //cout << "kmer = " << kmer << endl;
  string kmer_digit = Integer2KmerDigit(kmer);
  //cout << kmer_digit << endl;
  for (int i = 0; i < nNumNode; i++) {
    for (uint32_t j = 0; j < vLinkNodeSize[i]; j++) {
      if (vNodeLabel[vLinkList[i][j]] == 0) {
        vEdgeWeight[i][j] = 0;
        continue;
      }
      uint32_t id1 = vNodeLabel[vLinkList[i][j]] - 1;
      uint32_t id2 = kmer_digit[(i + ALPHABETSIZE - 1) / ALPHABETSIZE] - 48;
      vEdgeWeight[i][j] = REDUCEDBLOSUM62[id1][id2];
    }
  }
#ifdef debuggraph1
  //cout << "query=" << querySeq << endl;
  cout << "nNumNode = " << nNumNode << endl;
  for (int i = 0; i < nNumNode; i++) {
    for (uint32_t j = 0; j < vLinkNodeSize[i]; j++) {
      cout << "(" << i << ", " << vLinkList[i][j] << ", "
      << NODELABEL[vNodeLabel[vLinkList[i][j]]] << ", " << vEdgeWeight[i][j]
      << ")";
    }
    cout << endl;
  }
#endif
//////////////////////////////////////////////////////////
  for (int i = 0; i < nNumNode; i++) {
    adj[i].nxt = -1;
    invAdj[i].nxt = -1;
  }
  adjN = invAdjN = nNumNode;

  for (int i = 0; i < nNumNode; i++) {
    for (uint32_t j = 0; j < vLinkNodeSize[i]; j++) {
      int a, b;
      double c;
      a = i;
      b = vLinkList[i][j];
      c = vEdgeWeight[a][j];
      insert(adj, adjN, a, b, c);
      insert(invAdj, invAdjN, b, a, c);
    }
  }
  for (int i = 0; i < nNumNode; i++) {
    dp0[i] = 0;
    pre[i][0] = 0;
    while (!cand[i].empty())
      cand[i].pop();
  }
}

void KNearestNeighbor::BuildDAG() {
  /* construct nodes, there are totally 20*len + 2 nodes, the first node is source node,
   * and the last node is the destination node */
  int len = HASHLEN, nodeID = 0;
  /* set nodes label */
  vNodeLabel[nodeID++] = 0;
  for (int i = 0; i < len; i++) {
    for (int j = 0; j < ALPHABETSIZE; j++) {
      vNodeLabel[nodeID++] = j + 1;
    }
  }
  vNodeLabel[nodeID++] = 0;

  /* set edges */
  nNumNode = nodeID;
  int size = nodeID - ALPHABETSIZE - 1, base = 0;
  for (int i = 0; i < size; i++) {
    if (i % ALPHABETSIZE == 1)
      base += ALPHABETSIZE;
    for (int j = 1; j <= ALPHABETSIZE; j++) {
      vLinkList[i][j - 1] = base + j;
    }
    vLinkNodeSize[i] = ALPHABETSIZE;
  }

  for (int i = size; i < nNumNode - 1; i++) {
    vLinkList[i][0] = nNumNode - 1;
    vLinkNodeSize[i] = 1;
  }
  vLinkNodeSize[nNumNode - 1] = 0;

#ifdef TEST
  for(int i = 0;i < nNumNode;i++) {
    cout << " ii = " << i << endl;
    for(int j = 0;j < vLinkNodeSize[i];j++) {
      cout << vLinkList[i][j] << " ";
    }
    cout << endl;
  }
#endif

}

}  // namespace nearest_kmer