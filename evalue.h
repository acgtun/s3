#ifndef EVALUE_H_
#define EVALUE_H_

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
  Evalue(const uint64_t& _database_size,
         const uint32_t& _num_of_sequences_in_database)
      : database_size(_database_size),
        num_of_sequences_in_database(_num_of_sequences_in_database) {
    //cout << "database size = " << database_size << endl;
    //cout << "num of sequence = " << num_of_sequences_in_database << endl;

    Option::GetOption("-evalue", evalue_threshold, 10);
    log_database_size = log(database_size);
    logE = log(evalue_threshold);

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
  void UpdateValues(const uint32_t& query_len);

  double score_by_evalue[2];

 private:
  double evalue_threshold;
  double logE;

  uint64_t database_size;
  double log_database_size;
  uint32_t num_of_sequences_in_database;

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
