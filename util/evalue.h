#ifndef EVALUE_H_
#define EVALUE_H_

#include "bio_util.h"
#include "option.h"
#include <math.h>

const double lambda_gapped = 0.267;
const double K_gapped = 0.041;
const double H_gapped = 0.140;
const double logK_gapped = log(K_gapped);

const double lambda_ungapped = 0.318;
const double K_ungapped = 0.134;
const double H_ungapped = 0.401;
const double logK_ungapped = log(K_ungapped);

const double log_2 = log(2.0);

class Evalue {
 public:
  Evalue(const usint64_t& database_size,
         const usint32_t& num_of_sequences_in_database)
      : database_size_(database_size),
        num_of_sequences_in_database_(num_of_sequences_in_database) {
    cout << "database size = " << database_size_ << endl;
    cout << "num of sequence = " << num_of_sequences_in_database << endl;

    Option::GetOption("-evalue", evalue_threshold_, 10);
    log_database_size_ = log(database_size_);
    logE = log(evalue_threshold_);

    expected_HSP_length_ungapped = 1;
    expected_HSP_length_gapped = 1;

    effective_length_query_ungapped = 1;
    effective_length_database_ungapped = 1;

    effective_length_query_gapped = 1;
    effective_length_database_gapped = 1;

    KMN_expected_ungapped = 1;
    KMN_expected_gapped = 1;
  }
  ~Evalue() {
  }

  double GetBitScore(const int& score, const int& gap) const;
  double GetEvalue(const int& score, const int& gap) const;
  void UpdateValues(const usint32_t& query_len);

  double score_by_evalue_[2];

 private:
  double evalue_threshold_;
  double logE;

  usint64_t database_size_;
  double log_database_size_;
  usint32_t num_of_sequences_in_database_;

  double expected_HSP_length_ungapped;
  double expected_HSP_length_gapped;

  double effective_length_query_ungapped;
  double effective_length_database_ungapped;

  double effective_length_query_gapped;
  double effective_length_database_gapped;

  double KMN_expected_ungapped;
  double KMN_expected_gapped;
};

#endif /* EVALUE_H_ */
