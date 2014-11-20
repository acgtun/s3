#include "global_alignment.h"
#include "./../util/option.h"

#include <vector>
#include <limits.h>
#include "./../util/bio_util.h"

namespace global_alignment {
#define MAX_ALIGNMENT_LEN 6000

GlobalAlignment::GlobalAlignment() {
  Option::GetOption("-gopen", gapopen_, -11);
  Option::GetOption("-gext", gapextension_, -1);

  s.resize(MAX_ALIGNMENT_LEN);
  g.resize(MAX_ALIGNMENT_LEN);
  h.resize(MAX_ALIGNMENT_LEN);
  for (usint32_t i = 0; i < MAX_ALIGNMENT_LEN; i++) {
    s[i].resize(MAX_ALIGNMENT_LEN);
    g[i].resize(MAX_ALIGNMENT_LEN);
    h[i].resize(MAX_ALIGNMENT_LEN);
  }
}

GlobalAlignment::~GlobalAlignment() {

}

int GlobalAlignment::RunGlobalAlignment(const char * U, const char * V,
                                        const usint32_t & n,
                                        const usint32_t & m) {
#ifdef TESTCODE
  for (usint32_t i = 0; i < n; i++) {
    cout << U[i];
  }
  cout << endl;
  for (usint32_t i = 0; i < m; i++) {
    cout << V[i];
  }
  cout << endl;
  //start_t_ = clock();
#endif
  if (n >= MAX_ALIGNMENT_LEN || m >= MAX_ALIGNMENT_LEN) {
    ERROR_INFO("The length of the two sequences is too long.");
    return 0;
  }

  s[0][0] = 0;
  g[0][0] = INT_MIN + 10000;
  h[0][0] = INT_MIN + 10000;
  for (usint32_t i = 1; i <= m; i++) {
    s[0][i] = gapopen_ + i * gapextension_;
    g[0][i] = INT_MIN + 10000;
    h[0][i] = INT_MIN + 10000;
  }
  for (usint32_t i = 1; i <= n; i++) {
    s[i][0] = gapopen_ + i * gapextension_;
    g[i][0] = INT_MIN + 10000;
    h[i][0] = INT_MIN + 10000;
  }

  int stmp;
  for (usint32_t i = 1; i <= n; i++) {
    for (usint32_t j = 1; j <= m; j++) {
      g[i][j] = max(s[i][j - 1] + gapopen_ + gapextension_,
                    g[i][j - 1] + gapextension_);
      h[i][j] = max(s[i - 1][j] + gapopen_ + gapextension_,
                    h[i - 1][j] + gapextension_);

      stmp = s[i - 1][j - 1]
          + BLOSUM62[base[U[i - 1] - 'A']][base[V[j - 1] - 'A']];
      s[i][j] = max(max(g[i][j], h[i][j]), stmp);
      //cout << s[i][j] << " ";
    }
    //cout << endl;
  }
  //sum_time_ += (clock() - start_t_);

  // cout << "alignment score: " << s[n][m] << endl;
  return s[n][m];
}

}  // namespace local_alignment
