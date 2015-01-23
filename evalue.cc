#include "evalue.h"

double Evalue::GetBitScore(const int& score, const int& gap) const {
  if (gap == 0) {
    return (score * lambda_ungapped - logK_ungapped) / log_2;
  }
  return (score * lambda_gapped - logK_gapped) / log_2;
}

double Evalue::GetEvalue(const int& score, const int& gap) const {
  if (gap == 0) {
    return KMN_expected_ungapped * exp(-1.0 * lambda_ungapped * score);
    cout << KMN_expected_ungapped << " " << lambda_ungapped << " " << score
        << endl;
  }

  return KMN_expected_gapped * exp(-1.0 * lambda_gapped * score);
}

void Evalue::UpdateValues(const uint32_t& query_len) {
  expected_HSP_length_ungapped = (logK_ungapped + log_database_size
      + log(query_len)) / H_ungapped;
  expected_HSP_length_gapped =
      (logK_gapped + log_database_size + log(query_len)) / H_gapped;

  effective_length_query_ungapped = query_len - expected_HSP_length_ungapped;
  if (effective_length_query_ungapped < 1 / K_ungapped) {
    effective_length_query_ungapped = 1 / K_ungapped;
  }

  effective_length_query_gapped = query_len - expected_HSP_length_gapped;
  if (effective_length_query_gapped < 1 / K_gapped) {
    effective_length_query_gapped = 1 / K_gapped;
  }

  effective_length_database_ungapped = database_size
      - num_of_sequences_in_database * expected_HSP_length_ungapped;
  if (effective_length_database_ungapped < 1 / K_ungapped) {
    effective_length_database_ungapped = 1 / K_ungapped;
  }

  effective_length_database_gapped = database_size
      - num_of_sequences_in_database * expected_HSP_length_gapped;
  if (effective_length_database_gapped < 1 / K_gapped) {
    effective_length_database_gapped = 1 / K_gapped;
  }

  score_by_evalue[0] = (log(effective_length_database_ungapped)
      + log(effective_length_query_ungapped) + logK_ungapped - logE)
      / lambda_ungapped;
  score_by_evalue[1] = (log(effective_length_database_gapped)
      + log(effective_length_query_gapped) + logK_gapped - logE)
      / lambda_gapped;

  KMN_expected_ungapped = K_ungapped * effective_length_database_ungapped
      * effective_length_query_ungapped;

  KMN_expected_gapped = K_gapped * effective_length_database_gapped
      * effective_length_query_gapped;
}
