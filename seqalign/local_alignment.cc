#include "local_alignment.h"

#include <vector>
#include "./../util/bio_util.h"

namespace local_alignment {

LocalAlignment::LocalAlignment(const Evalue* evalue, const usint32_t& max_rows,
                               const usint32_t& max_cols)
    : evalue_(evalue),
      max_rows_(max_rows),
      max_cols_(max_cols) {

  Option::GetOption("-gopen", gapopen_, -11);
  Option::GetOption("-gext", gapextension_, -1);

  s.resize(max_rows + 3);
  l.resize(max_rows + 3);

  g.resize(max_rows + 3);
  h.resize(max_rows + 3);

  for (usint32_t i = 0; i <= max_rows; i++) {
    s[i].resize(max_cols + 3);
    l[i].resize(max_cols + 3);

    g[i].resize(max_cols + 3);
    h[i].resize(max_cols + 3);
  }
  MEMORY_ALLOCATE_CHECK(
      rU = (char *) malloc(sizeof(char) * (max_rows + max_cols + 2)));
  MEMORY_ALLOCATE_CHECK(
      rV = (char *) malloc(sizeof(char) * (max_rows + max_cols + 2)));
  MEMORY_ALLOCATE_CHECK(
      midline = (char *) malloc(sizeof(char) * (max_rows + max_cols + 2)));

  alignScore = 0;
  nIdentity = 0;

  m = 0;
  n = 0;
  rs = 0;
  start_t_ = 0;
  sum_time_ = 0;
}

LocalAlignment::~LocalAlignment() {
  INFO("Release Local Alignment memory...");
  free(rU);
  free(rV);
  free(midline);

  printf("--INFO-- Local Alignment takes %.3lf seconds.\n",
         (double) sum_time_ / CLOCKS_PER_SEC);
}

void LocalAlignment::stringReverse(char * str, const int & n) {
  char c;
  for (int i = 0; i < n / 2; i++) {
    c = str[i];
    str[i] = str[n - i - 1];
    str[n - i - 1] = c;
  }
}

int LocalAlignment::MaxOfFour(const int & s1, const int & s2, const int & s3,
                              const int & s4) {
  /*if two of them are equal, then there are more than one optimal path*/
  return max(max(max(s1, s2), s3), s4);
}

char LocalAlignment::Direction(const int & s1, const int & s2, const int & s3,
                               const int & s4) {
  if (s4 == s2)
    return UP;
  else if (s4 == s3)
    return LEFT;
  else if (s4 == s1)
    return DIAG;
  else
    //(s4 == 0)
    return STOPorSTART;
}

void LocalAlignment::DisplayAlignment() {
  nIdentity = 0;
  for (int i = 0; i < rs; i++) {
    if (rU[i] == rV[i]) {
      nIdentity++;
    }
  }

  for (int i = 0; i < rs; i++) {
    if (rU[i] == rV[i]) {
      midline[i] = rV[i];
    } else if (rU[i] != '-' && rV[i] != '-'
        && BLOSUM62[base[rU[i] - 'A']][base[rV[i] - 'A']] > 0) {
      midline[i] = '+';
    } else {
      midline[i] = ' ';
    }
  }
  midline[rs] = 0;

  //cout << rU << endl;
  //cout << midline << endl;
 // cout << rV << endl;
  //cout << "Alignment score: " << alignScore << endl;
 // cout << "Identity: " << (double) nIdentity / rs * 100 << "%" << endl;
  //cout << "evalue: " << e_value << endl;
}

bool LocalAlignment::RunLocalAlignment(const char * U, const char * V,
                                       const usint32_t & Ul,
                                       const usint32_t & Vl, M8Results& res) {
  start_t_ = clock();
  if (Ul > max_rows_ || Vl > max_cols_) {
    printf("The length of the two sequences is too long.");
    return false;
  }
  n = Ul;
  m = Vl;
  s[0][0] = 0;
  g[0][0] = INT_MIN + 10000;
  h[0][0] = INT_MIN + 10000;

  for (usint32_t i = 1; i <= m; i++) {
    s[0][i] = 0;
    l[0][i] = LEFT;
    g[0][i] = INT_MIN + 10000;
    h[0][i] = INT_MIN + 10000;
  }
  for (usint32_t i = 1; i <= n; i++) {
    s[i][0] = 0;
    l[0][i] = UP;
    g[i][0] = INT_MIN + 10000;
    h[i][0] = INT_MIN + 10000;
  }
  int stmp, s4;
  int max_score = 0;
  usint32_t max_i = 0, max_j = 0;
  for (usint32_t i = 1; i <= n; i++) {
    for (usint32_t j = 1; j <= m; j++) {
      g[i][j] = max(s[i][j - 1] + gapopen_ + gapextension_,
                    g[i][j - 1] + gapextension_);
      h[i][j] = max(s[i - 1][j] + gapopen_ + gapextension_,
                    h[i - 1][j] + gapextension_);
      stmp = s[i - 1][j - 1]
          + BLOSUM62[base[U[i - 1] - 'A']][base[V[j - 1] - 'A']];
      s4 = MaxOfFour(0, stmp, h[i][j], g[i][j]);
      s[i][j] = s4;
      l[i][j] = Direction(stmp, h[i][j], g[i][j], s4);
      if (s[i][j] > max_score) {
        max_score = s[i][j];
        max_i = i;
        max_j = j;
      }
    }
  }

  int p = max_i, q = max_j;
  rs = 0;
  int gap = 0;
  //cout << "--------------------" << endl;
  while (p > 0 && q > 0) {  // trace back from s[n][m] to s[0][0]
    if (l[p][q] == DIAG) {
      rU[rs] = U[p - 1];
      rV[rs] = V[q - 1];
      rs++;
      p = p - 1;
      q = q - 1;
    } else if (l[p][q] == UP) {
      rU[rs] = U[p - 1];
      rV[rs] = '-';
      rs++;
      p = p - 1;
      gap = 1;
    } else if (l[p][q] == LEFT) {
      rU[rs] = '-';
      rV[rs] = V[q - 1];
      rs++;
      q = q - 1;
      gap = 1;
    } else {
      break;
    }
  }
  rU[rs] = 0;
  rV[rs] = 0;
  stringReverse(rU, rs);
  stringReverse(rV, rs);

  alignScore = max_score;
  //DisplayAlignment();
  if (alignScore >= evalue_->score_by_evalue_[gap]) {
    e_value = evalue_->GetEvalue(alignScore, gap);
    DisplayAlignment();
    res = M8Results("", nIdentity, rs, 0, 0, p + 1, max_i, q + 1, max_j,
                    e_value, evalue_->GetBitScore(alignScore, gap));
    return true;
  }
  sum_time_ += (clock() - start_t_);
  return false;
}

}  // namespace local_alignment
