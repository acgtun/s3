#ifndef HITS_MAXTRIX_H_
#define HITS_MAXTRIX_H_

#include "maximal_path.h"
#include "./util/option.h"
#include "./index/index.h"
#include "./util/fasta_file.h"
#include <stack>
#include <list>
#include "./util/bio_util.h"
#include "./util/evalue.h"
#include "./seqalign/local_alignment.h"

using maximal_path::CSegHit;
using maximal_path::MaximalPath;
using local_alignment::LocalAlignment;

namespace hits_matrix {

class HitsMatrix {
 public:
  HitsMatrix(const Evalue* evalue, const usint32_t& max_query_length,
             const usint32_t& max_protein_length)
      : evalue_(evalue),
        maximal_path_(evalue),
        local_alignment_(evalue, max_query_length, max_protein_length) {
    Option::GetOption("-D", D, 5);
    Option::GetOption("-diagdiff", diags_diff, 50);
    Option::ChkStrExist("-local", using_local_alignment);

    query_ = NULL;
    protein_ = NULL;
  }
  ~HitsMatrix() {
  }

  void AnalyzeMatrix(const char* query, const char* protein,
                     const string& protein_name, const vector<CSegHit>& hits,
                     vector<M8Results>& aligned_results);

 private:
  void ChainFilter(map<int, vector<CSegHit> >& seg,
                   map<int, vector<CSegHit> >& seg_new);
  void SelectDiags(map<int, vector<CSegHit> >& seg,
                   vector<M8Results>& aligned_results);
  void MatlabShowDiags(map<int, vector<CSegHit> >& seg);
  int DiagID(const usint32_t & qs, const usint32_t & ps);
  void UngappedExtension(map<int, vector<CSegHit> >& seg);
  void SplitToDiag(const vector<CSegHit>& hits,
                   map<int, vector<CSegHit> >& seg);
  int SegScore(const CSegHit& v);
  void ShowHits(map<int, vector<CSegHit> >& seg);

  const char* query_;
  const char* protein_;
  string protein_name_;

  const Evalue* evalue_;
  MaximalPath maximal_path_;
  LocalAlignment local_alignment_;

  int D;
  int diags_diff;
  bool using_local_alignment;
};

}
#endif /* HITS_MAXTRIX_H_ */
