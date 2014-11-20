#include "fasta_file.h"

void FastaFile::SetFilePath(const string& file_path) {
  file_path_ = file_path;
  ReadFastaFile();
}

string FastaFile::GetFilePath() const {
  return file_path_;
}

void FastaFile::ReadFastaFile() {
  INFO("Read Fasta file:", file_path_.c_str());
  FILE* fin = fopen(file_path_.c_str(), "rb");
  FILE_OPEN_CHECK(fin);
  fseek(fin, 0, SEEK_END);
  file_size_ = ftell(fin);
  INFO("file_size", file_size_);
  MEMORY_ALLOCATE_CHECK(
      file_string_ = (char* ) malloc(sizeof(char) * (file_size_ + 1)));
  fseek(fin, 0, SEEK_SET);
  FREAD_CHECK(fread(file_string_, 1, file_size_, fin), file_size_);
  fclose(fin);

  sequences_ = NULL;
  num_of_sequences_ = 0;
  num_of_characters_ = 0;

  GetNumOfSequences();
  AnalyzeFile();
}

FastaFile::FastaFile(const string& file_path) {
  file_path_ = file_path;
  file_string_ = NULL;
  file_size_ = 0;

  sequences_ = NULL;
  num_of_sequences_ = 0;
  num_of_characters_ = 0;
  max_sequence_length_ = 100000;
  ReadFastaFile();
}

FastaFile::FastaFile() {
  file_path_ = "";
  file_string_ = NULL;
  file_size_ = 0;

  sequences_ = NULL;
  num_of_sequences_ = 0;
  num_of_characters_ = 0;
  max_sequence_length_ = 100000;
}

FastaFile::~FastaFile() {
  for (usint32_t i = 0; i < num_of_sequences_; i++) {
    free(sequences_[i]);
  }
  free(sequences_);
  sequences_names_.clear();
}

usint64_t FastaFile::GetLineFromString(const char* strVal, char* strRet) {
  usint64_t i;
  bool tag = 0;
  usint64_t j = 0;
  for (i = 0; strVal[i] != 0; i++) {
    if (strVal[i] == ' ')
      tag = 1;
    if (0xA == strVal[i] || 0xD == strVal[i]) {
      break;
    }
    if (tag == 0) {
      strRet[j] = strVal[i];
      j++;
    }
  }

  strRet[j] = 0;
  return i;
}

void FastaFile::GetNumOfSequences() {
  char strRet[MAX_LINE_LEN];
  for (usint64_t i = 0; i < file_size_; i++) {
    if (file_string_[i] == '>') {
      num_of_sequences_++;
    }
    i += GetLineFromString(&file_string_[i], strRet);
  }
}

void FastaFile::AnalyzeFile() {
  MEMORY_ALLOCATE_CHECK(
      sequences_ = (char** ) malloc(sizeof(char *) * num_of_sequences_));
  sequences_names_.resize(num_of_sequences_);

  uint32_t id = 0;
  string str_seq_;
  num_of_characters_ = 0;
  max_sequence_length_ = 0;
  str_seq_.clear();
  char strRet[MAX_LINE_LEN];
  for (usint64_t i = 0; i < file_size_; i++) {
    if (file_string_[i] == '>') {
      if (str_seq_.size() != 0) {
        MEMORY_ALLOCATE_CHECK(
            sequences_[id] = (char * ) malloc(
                sizeof(char) * (str_seq_.size() + 1)));
        strcpy(sequences_[id], str_seq_.c_str());
        num_of_characters_ += str_seq_.size();
        if (str_seq_.size() > max_sequence_length_) {
          max_sequence_length_ = str_seq_.size();
        }
        str_seq_.clear();
        id++;
      }
      i += GetLineFromString(&file_string_[i], strRet);
      sequences_names_[id] = string(strRet).substr(1);
    } else if (AA20.find_first_of(file_string_[i]) != string::npos) {
      str_seq_ += toupper(file_string_[i]);
    } else if (isalpha(file_string_[i])) {
      str_seq_ += AA20[rand() % 20];  //todo
    }
  }
  if (str_seq_.size() != 0) {
    MEMORY_ALLOCATE_CHECK(
        sequences_[id] = (char * ) malloc(
            sizeof(char) * (str_seq_.size() + 1)));
    strcpy(sequences_[id], str_seq_.c_str());
    num_of_characters_ += str_seq_.size();
    if (str_seq_.size() > max_sequence_length_) {
      max_sequence_length_ = str_seq_.size();
    }
    str_seq_.clear();
    id++;
  }
  if (id != num_of_sequences_) {
    ERROR_INFO("The number of Sequences has error~!");
  }

  free(file_string_);

  //////////////////////////////////////////////////////
  FILE * fout = fopen("proteins.txt", "w");
  for (usint32_t i = 0; i < num_of_sequences_; i++) {
    fprintf(fout, "%d:%s\n", i, sequences_names_[i].c_str());
    fprintf(fout, "%s\n", sequences_[i]);
  }
  fclose(fout);
  ///////////////////////////////////////////////////////
}
