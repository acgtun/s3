#include "buildindex.h"
#include "iofile.h"
#include "KlongestPathDP.h"
#include <set>
#include <stdlib.h>
#include <string>
#include <stdio.h>
#include <string.h>

#define generate_test_data
#define proDBSee

Index::Index(const Option& opt) :
    opt_(opt) {
      score_threshold_ = 10;
    }

Index::~Index() {
  for (usint32_t i = 0; i < proDB_.nNumPro; i++) {
    free(proDB_.vProSeq[i]);
  }
  free(proDB_.vProSeq);
  proDB_.vProName.clear();
  proDB_.vProSize.clear();

  free(kmer_nearest_.counter);
  free(kmer_nearest_.index);

  free(kmer_position_.counter);
  free(kmer_position_.index_proID);
  free(kmer_position_.index_proPos);
}

void Index::GetThresholdScore() {
  INFO("Get Threshold Score...");
  map<int, usint64_t> org_count;
  for (usint32_t i = 0; i < AA20.size(); i++) {
    for (usint32_t j = 0; j < AA20.size(); j++) {
      org_count[BLOSUM62[i][j]]++;
    }
  }

  map<int, usint64_t> new_count;
  map<int, usint64_t> old_count(org_count);
  for (int i = 1; i < HASHAALEN; i++) {
    new_count.clear();
    for (map<int, usint64_t>::iterator it1 = old_count.begin(); it1 != old_count.end(); it1++) {
      for (map<int, usint64_t>::iterator it2 = org_count.begin(); it2 != org_count.end(); it2++) {
        new_count[it1->first + it2->first] += it1->second * it2->second;
      }
    }
    old_count.clear();
    old_count = new_count;
  }

  usint64_t sum = 0;
  INFO("Calculate Threshold Score...");
  for (map<int, usint64_t>::iterator it = old_count.begin(); it != old_count.end(); it++) {
    //cout << it->first << " " << it->second << endl;
    sum += it->second;
  }
  //cout << "sum = " << sum << endl;

  ////////////////////////////////////////////////////////////////////
  cout << proDB_.sizeAA << endl;
  double p = (double) opt_.expected_hits / (double) proDB_.sizeAA;
  INFO("The probability for each pair that has score larger than T is ", p);
  usint64_t threshold = 0;
  for (map<int, usint64_t>::reverse_iterator it = old_count.rbegin(); it != old_count.rend(); it++) {
    cout << it->first << " " << it->second << endl;
    threshold += it->second;
    if ((double) threshold / (double) sum > p) {
      score_threshold_ = it->first;
      break;
    }
  }
  INFO("The threshold score is", score_threshold_);
}

void Index::CountNearestKmer() {
  INFO("Count Nearest Kmer...");
  kmer_nearest_.nSizeCounter = GetPower(20, HASHAALEN);
  INFO("The size of Kmer Nearest Counter is ", kmer_nearest_.nSizeCounter);
  MEMORY_ALLOCATE_CHECK(
      kmer_nearest_.counter = (usint64_t * ) malloc(sizeof(usint64_t) * (kmer_nearest_.nSizeCounter + 1)));
  memset(kmer_nearest_.counter, 0x00, (kmer_nearest_.nSizeCounter + 1));

  usint32_t hashValueCand = 0;
  char querySeq[HASHAALEN + 1], strCandPep[HASHAALEN + 1];
  for (usint32_t i = 0; i < kmer_position_.nSizeCounter; i++) {
    if (kmer_position_.counter[i] == kmer_position_.counter[i + 1])
      continue;
    DeCodeAA(i, querySeq);
    m_kLongestPath_.UpdateWeight(querySeq);
    bool bPath = true;
    int pathID = 0, rawScore = 0;
    while (bPath && pathID < MAX_K) {
      bPath = m_kLongestPath_.FindCandidiate(strCandPep, pathID, rawScore);
      if (pathID == 0 && rawScore < score_threshold_) {
        if (strcmp(strCandPep, querySeq) != 0) {
          printf("CAND=%s\n", strCandPep);
          printf("QUERY=%s\n", querySeq);
          cout << "i = " << i << " " << "rawScore = " << rawScore << endl;
          ERROR_INFO("The Find k Longest Path has error");
        }
        kmer_nearest_.counter[i]++;
        break;
      }
      if (rawScore < score_threshold_)
        break;
      hashValueCand = GetHashValue(strCandPep);
      kmer_nearest_.counter[hashValueCand]++;
      pathID++;
    }
  }

  // 2^32 / 20^6
  // Here we should be careful about the sum since the sum may more
  // than 2^32 if the average nubmer of nearest for each query is
  // more than certain number.
  for (usint32_t i = 1; i <= kmer_nearest_.nSizeCounter; i++) {
    kmer_nearest_.counter[i] += kmer_nearest_.counter[i - 1];
  }
  kmer_nearest_.nSizeIndex = kmer_nearest_.counter[kmer_nearest_.nSizeCounter];
  INFO("The size of Kmer Nearest Index is ", kmer_nearest_.nSizeIndex);

  for (usint32_t i = kmer_nearest_.nSizeCounter - 1; i >= 1; i--) {
    kmer_nearest_.counter[i] = kmer_nearest_.counter[i - 1];
  }
  kmer_nearest_.counter[0] = 0;
}

void Index::BuildNearestKmer() {
  INFO("Build Nearest Kmer...");
  MEMORY_ALLOCATE_CHECK(kmer_nearest_.index = (usint32_t * ) malloc(sizeof(usint32_t) * kmer_nearest_.nSizeIndex));
  usint32_t hashValueCand = 0;
  char querySeq[HASHAALEN + 1], strCandPep[HASHAALEN + 1];
  for (usint32_t i = 0; i < kmer_position_.nSizeCounter; i++) {
    if (kmer_position_.counter[i] == kmer_position_.counter[i + 1])
      continue;
    DeCodeAA(i, querySeq);
    m_kLongestPath_.UpdateWeight(querySeq);
    bool bPath = true;
    int pathID = 0, rawScore = 0;
    while (bPath && pathID < MAX_K) {
      bPath = m_kLongestPath_.FindCandidiate(strCandPep, pathID, rawScore);
      if (pathID == 0 && rawScore < score_threshold_) {
        if (strcmp(strCandPep, querySeq) != 0) {
          ERROR_INFO("The Find k Longest Path has error");
        }
        kmer_nearest_.index[kmer_nearest_.counter[i]++] = i;
        break;
      }
      if (rawScore < score_threshold_)
        break;
      hashValueCand = GetHashValue(strCandPep);
      kmer_nearest_.index[kmer_nearest_.counter[hashValueCand]++] = i;
      pathID++;
    }
  }

  for (usint32_t i = kmer_nearest_.nSizeCounter - 1; i >= 1; i--) {
    kmer_nearest_.counter[i] = kmer_nearest_.counter[i - 1];
  }
  kmer_nearest_.counter[0] = 0;
}

void Index::CountKmerPosition() {
  INFO("Count Kmer Position...");
  kmer_position_.nSizeCounter = GetPower(20, HASHAALEN);
  INFO("The size of Kmer Position Counter is", kmer_position_.nSizeCounter);
  MEMORY_ALLOCATE_CHECK(
      kmer_position_.counter = (usint64_t * ) malloc(sizeof(usint64_t) * (kmer_position_.nSizeCounter + 1)));
  memset(kmer_position_.counter, 0x00, (kmer_position_.nSizeCounter + 1));

  usint32_t size = 0, hashValue = 0;
  for (usint32_t i = 0; i < proDB_.nNumPro; i++) {
    if (proDB_.vProSize[i] < HASHAALEN) {
      INFO("The length of the protein is less than 6, ignore it.", proDB_.vProName[i]);
      continue;
    }
    size = proDB_.vProSize[i] - HASHAALEN;
    for (usint32_t j = 0; j <= size; j++) {
      hashValue = GetHashValue(&(proDB_.vProSeq[i][j]));
      kmer_position_.counter[hashValue]++;
    }
  }

  for (usint32_t i = 1; i <= kmer_position_.nSizeCounter; i++) {
    kmer_position_.counter[i] += kmer_position_.counter[i - 1];
  }
  kmer_position_.nSizeIndex = kmer_position_.counter[kmer_position_.nSizeCounter];
  INFO("The size of Kmer Position Index is", kmer_position_.nSizeIndex);

  for (usint32_t i = kmer_position_.nSizeCounter - 1; i >= 1; i--) {
    kmer_position_.counter[i] = kmer_position_.counter[i - 1];
  }
  kmer_position_.counter[0] = 0;
}

void Index::BuildKmerPosition() {
  INFO("Build Kmer Position...");
  MEMORY_ALLOCATE_CHECK(
      kmer_position_.index_proID = (usint32_t * ) malloc(sizeof(usint32_t) * kmer_position_.nSizeIndex));
  MEMORY_ALLOCATE_CHECK(
      kmer_position_.index_proPos = (usint32_t * ) malloc(sizeof(usint32_t) * kmer_position_.nSizeIndex));
  usint32_t size = 0, hashValue = 0;
  for (usint32_t i = 0; i < proDB_.nNumPro; i++) {
    if (proDB_.vProSize[i] < HASHAALEN) {
      INFO("The length of the protein is less than 6, ignore it.", proDB_.vProName[i]);
      continue;
    }
    size = proDB_.vProSize[i] - HASHAALEN;
    for (usint32_t j = 0; j <= size; j++) {
      hashValue = GetHashValue(&(proDB_.vProSeq[i][j]));
      //test////////////////////////////////////////////
      //char testchar[HASHAALEN + 1];
      //DeCodeAA(hashValue, testchar);
      //cout << hashValue << " " << proDB_.vProSeq[i][j] << proDB_.vProSeq[i][j + 1] << proDB_.vProSeq[i][j + 2]
      //<< proDB_.vProSeq[i][j + 3] << proDB_.vProSeq[i][j + 4] << proDB_.vProSeq[i][j + 5] << " "
      //<< testchar << endl;
      //test////////////////////////////////////////////
      kmer_position_.index_proID[kmer_position_.counter[hashValue]] = i;
      kmer_position_.index_proPos[kmer_position_.counter[hashValue]] = j;
      kmer_position_.counter[hashValue]++;
    }
  }

  for (usint32_t i = kmer_position_.nSizeCounter - 1; i >= 1; i--) {
    kmer_position_.counter[i] = kmer_position_.counter[i - 1];
  }
  kmer_position_.counter[0] = 0;
}

void Index::AnalyzeProteinDataBase(const char * strRef, const usint64_t & refLen) {
  INFO("Analyze the Protein Database...");
  proDB_.nNumPro = 0;
  char strRet[MAX_LINE_LEN]; // proName[MAX_LINE_LEN];
  for (usint64_t i = 0; i < refLen; i++) {
    if (strRef[i] == '>') {
      proDB_.nNumPro++;
    }
    i += GetLineFromString(&strRef[i], strRet);
  }

  MEMORY_ALLOCATE_CHECK(proDB_.vProSeq = (char ** ) malloc(sizeof(char *) * proDB_.nNumPro));
  proDB_.vProName.resize(proDB_.nNumPro);
  proDB_.vProSize.resize(proDB_.nNumPro);

  uint32_t nProID = 0;
  string strProSeq;
  proDB_.sizeAA = 0;
  strProSeq.clear();
  srand(time(NULL));
  for (usint64_t i = 0; i < refLen; i++) {
    if (strRef[i] == '>') {
      if (strProSeq.size() != 0) {
        MEMORY_ALLOCATE_CHECK(proDB_.vProSeq[nProID] = (char * ) malloc(sizeof(char) * (strProSeq.size() + 1)));
        strcpy(proDB_.vProSeq[nProID], strProSeq.c_str());
        proDB_.vProSize[nProID] = strProSeq.size();
        proDB_.sizeAA += proDB_.vProSize[nProID];
        strProSeq.clear();
        nProID++;
      }
      else {
        if(nProID > 0) {
          INFO("Something error here~!");
        }
      }
      //next protein
      i += GetLineFromString(&strRef[i], strRet);
      proDB_.vProName[nProID] = strRet;
    } else if (AA20.find_first_of(strRef[i]) != string::npos) {
      strProSeq += toupper(strRef[i]);
    } else if(isalpha(strRef[i])) {
      strProSeq += AA20[rand() % 20];
    }
  }
  if (strProSeq.size() != 0) {
    MEMORY_ALLOCATE_CHECK(proDB_.vProSeq[nProID] = (char * ) malloc(sizeof(char) * (strProSeq.size() + 1)));
    strcpy(proDB_.vProSeq[nProID], strProSeq.c_str());
    proDB_.vProSize[nProID] = strProSeq.size();
    proDB_.sizeAA += proDB_.vProSize[nProID];
    strProSeq.clear();
    nProID++;
  }
  if (nProID != proDB_.nNumPro) {
    cout << "nProID = " << nProID << endl;
    cout << "proDB_.nNumPro = " << proDB_.nNumPro << endl;
    ERROR_INFO("The number of Protein is error");
  }
  INFO("There are", proDB_.nNumPro, "proteins in the database");
  INFO("The total length of these proteins is", proDB_.sizeAA);
#ifdef proDBSee
  FILE * fout = fopen("proteins.txt", "w");
  for (usint32_t i = 0; i < proDB_.nNumPro; i++) {
    fprintf(fout, "%s\n", proDB_.vProName[i].c_str());
    fprintf(fout, "%s\n", proDB_.vProSeq[i]);
  }
  fclose(fout);

  ofstream ftd("test_database.fasta");
  for (usint32_t i = 0; i < 100; i++) {
    ftd << proDB_.vProName[i] << endl;
    ftd << proDB_.vProSeq[i] << endl;
  }
  ftd.close();
#endif
#ifdef generate_test_data
  srand(time(NULL));
  int testID = 0;
  ofstream fall("test_data/test_10000_all.fasta");
  for(int d = 0; d < 1000;d++) {
    char file_name[100];
    sprintf(file_name, "test_data/test_10000_%d.fasta", d);
    ftd.open(file_name);
    for (usint32_t i = 0; i < 10; i++) {
      usint32_t id, start_pos, len;
      while(1) {
        id = rand() % proDB_.nNumPro;
        start_pos = rand() % 1000;
        len = rand() % 1000 + 6;
        if (start_pos + len < proDB_.vProSize[id])
          break;
      }
      ftd << ">testprotein" << testID << " " << proDB_.vProName[id] << " " << start_pos << " " << len << endl;
      fall << ">testprotein" << testID << " " << proDB_.vProName[id] << " " << start_pos << " " << len << endl;
      testID++;
      //ftd << "size =  " <<  proDB_.vProSize[id] << endl;
      //ftd << "seq=" << proDB_.vProSeq[id] << endl;
      // gap
      for (usint32_t k = start_pos, l = 0; l < len; k++, l++) {
        if(l % 80 == 0 && l != 0) {
          ftd << endl;
          fall << endl;
        }
        int r = rand() % 250;
        if(r == 0) continue;
        ftd << proDB_.vProSeq[id][k];
        fall << proDB_.vProSeq[id][k];
      }
      ftd << endl;
      fall << endl;
    }
    ftd.close();
  }
  fall.close();
  exit(EXIT_SUCCESS);
#endif
}

void Index::BuildIndex() {
  char * strRef;
  INFO("Read protein database from", opt_.refFile);
  usint64_t refLen = ReadWholeFile(opt_.refFile, &strRef);
  AnalyzeProteinDataBase(strRef, refLen);
  free(strRef);
  GetThresholdScore();

  TIME_INFO(CountKmerPosition(), "Count Kmer Position");
  TIME_INFO(BuildKmerPosition(), "Build Kmer Position");

  TIME_INFO(CountNearestKmer(), "Count Neareast Kmer");
  TIME_INFO(BuildNearestKmer(), "Build Nearest Kmer");

  //TIME_INFO(TestIndex(), "Write Test Index");
  TIME_INFO(WriteIndex(), "Write Index");
}

void Index::WriteIndex() {
  FILE * fout = fopen(opt_.indexFile.c_str(), "wb");
  INFO("Write index to", opt_.indexFile);

  fwrite(&(proDB_.sizeAA), sizeof(usint64_t), 1, fout);
  fwrite(&(proDB_.nNumPro), sizeof(usint32_t), 1, fout);
  for (usint32_t i = 0; i < proDB_.nNumPro; i++) {
    usint32_t l = proDB_.vProName[i].size();
    fwrite(&l, sizeof(usint32_t), 1, fout);
    fwrite(&(proDB_.vProName[i][0]), sizeof(char), l, fout);
    fwrite(&(proDB_.vProSize[i]), sizeof(usint32_t), 1, fout);
    fwrite(proDB_.vProSeq[i], sizeof(char), proDB_.vProSize[i], fout);
  }

  fwrite(&(kmer_nearest_.nSizeCounter), sizeof(usint32_t), 1, fout);
  fwrite(&(kmer_nearest_.nSizeIndex), sizeof(usint64_t), 1, fout);
  fwrite(kmer_nearest_.counter, sizeof(usint64_t), kmer_nearest_.nSizeCounter + 1, fout);
  fwrite(kmer_nearest_.index, sizeof(usint32_t), kmer_nearest_.nSizeIndex, fout);

  fwrite(&(kmer_position_.nSizeCounter), sizeof(usint32_t), 1, fout);
  fwrite(&(kmer_position_.nSizeIndex), sizeof(usint64_t), 1, fout);
  fwrite(kmer_position_.counter, sizeof(usint64_t), kmer_position_.nSizeCounter + 1, fout);
  fwrite(kmer_position_.index_proID, sizeof(usint32_t), kmer_position_.nSizeIndex, fout);
  fwrite(kmer_position_.index_proPos, sizeof(usint32_t), kmer_position_.nSizeIndex, fout);

  fclose(fout);
}

void Index::TestIndex() {
  INFO("Write Test File...");
  ofstream fp("table_position.txt");
  char querySeq[HASHAALEN + 1];
  for (usint32_t i = 0; i < kmer_position_.nSizeCounter; i++) {
    if (kmer_position_.counter[i] == kmer_position_.counter[i + 1])
      continue;
    DeCodeAA(i, querySeq);
    fp << i << " " << querySeq << " " << kmer_position_.counter[i] << endl;
    for (usint32_t j = kmer_position_.counter[i]; j < kmer_position_.counter[i + 1]; j++) {
      fp << proDB_.vProName[kmer_position_.index_proID[j]] << " " << kmer_position_.index_proPos[j] << " ";
      fp << proDB_.vProSeq[kmer_position_.index_proID[j]][kmer_position_.index_proPos[j]]
          << proDB_.vProSeq[kmer_position_.index_proID[j]][kmer_position_.index_proPos[j] + 1]
          << proDB_.vProSeq[kmer_position_.index_proID[j]][kmer_position_.index_proPos[j] + 2]
          << proDB_.vProSeq[kmer_position_.index_proID[j]][kmer_position_.index_proPos[j] + 3]
          << proDB_.vProSeq[kmer_position_.index_proID[j]][kmer_position_.index_proPos[j] + 4]
          << proDB_.vProSeq[kmer_position_.index_proID[j]][kmer_position_.index_proPos[j] + 5] << endl;
    }
  }
  fp.close();

  ofstream fn("table_nearest.txt");
  char queryCand[HASHAALEN + 1];
  for (usint32_t i = 0; i < kmer_nearest_.nSizeCounter; i++) {
    if (kmer_nearest_.counter[i] == kmer_nearest_.counter[i + 1])
      continue;
    DeCodeAA(i, querySeq);
    fn << querySeq << " (";
    for (usint32_t j = kmer_nearest_.counter[i]; j < kmer_nearest_.counter[i + 1]; j++) {
      DeCodeAA(kmer_nearest_.index[j], queryCand);
      fn << queryCand << " ,";
    }
    fn << ")" << endl;
  }
  fn.close();

}
void Index::ReadIndex() {
  INFO("Read index from", opt_.refFile);
  FILE * fin = fopen(opt_.refFile.c_str(), "rb");
  FILE_OPEN_CHECK(fin);

  FREAD_CHECK(fread(&(proDB_.sizeAA), sizeof(usint64_t), 1, fin), 1);
  FREAD_CHECK(fread(&(proDB_.nNumPro), sizeof(usint32_t), 1, fin), 1);
  char proName[MAX_LINE_LEN];
  usint32_t l;
  MEMORY_ALLOCATE_CHECK(proDB_.vProSeq = (char ** ) malloc(sizeof(char *) * proDB_.nNumPro));
  for (usint32_t i = 0; i < proDB_.nNumPro; i++) {
    FREAD_CHECK(fread(&l, sizeof(usint32_t), 1, fin), 1);
    FREAD_CHECK(fread(proName, sizeof(char), l, fin), l);
    proName[l] = 0;
    proDB_.vProName.push_back(proName);
    FREAD_CHECK(fread(&l, sizeof(usint32_t), 1, fin), 1);
    proDB_.vProSize.push_back(l);
    MEMORY_ALLOCATE_CHECK(proDB_.vProSeq[i] = (char * ) malloc(sizeof(char) * (l + 1)));
    FREAD_CHECK(fread(&(proDB_.vProSeq[i][0]), sizeof(char), l, fin), l);
    proDB_.vProSeq[i][l] = 0;
  }

  // kmer_nearest_
  FREAD_CHECK(fread(&(kmer_nearest_.nSizeCounter), sizeof(usint32_t), 1, fin), 1);
  FREAD_CHECK(fread(&(kmer_nearest_.nSizeIndex), sizeof(usint64_t), 1, fin), 1);
  MEMORY_ALLOCATE_CHECK(
      kmer_nearest_.counter = (usint64_t * ) malloc(sizeof(usint64_t) * (kmer_nearest_.nSizeCounter + 1)));
  MEMORY_ALLOCATE_CHECK(kmer_nearest_.index = (usint32_t * ) malloc(sizeof(usint32_t) * kmer_nearest_.nSizeIndex));
  FREAD_CHECK(fread(kmer_nearest_.counter, sizeof(usint64_t), kmer_nearest_.nSizeCounter + 1, fin),
              kmer_nearest_.nSizeCounter + 1);
  FREAD_CHECK(fread(kmer_nearest_.index, sizeof(usint32_t), kmer_nearest_.nSizeIndex, fin), kmer_nearest_.nSizeIndex);

  // kmer_position_
  FREAD_CHECK(fread(&(kmer_position_.nSizeCounter), sizeof(usint32_t), 1, fin), 1);
  FREAD_CHECK(fread(&(kmer_position_.nSizeIndex), sizeof(usint64_t), 1, fin), 1);
  MEMORY_ALLOCATE_CHECK(
      kmer_position_.counter = (usint64_t * ) malloc(sizeof(usint64_t) * (kmer_position_.nSizeCounter + 1)));
  MEMORY_ALLOCATE_CHECK(
      kmer_position_.index_proID = (usint32_t * ) malloc(sizeof(usint32_t) * kmer_position_.nSizeIndex));
  MEMORY_ALLOCATE_CHECK(
      kmer_position_.index_proPos = (usint32_t * ) malloc(sizeof(usint32_t) * kmer_position_.nSizeIndex));
  FREAD_CHECK(fread(kmer_position_.counter, sizeof(usint64_t), kmer_position_.nSizeCounter + 1, fin),
              kmer_position_.nSizeCounter + 1);
  FREAD_CHECK(fread(kmer_position_.index_proID, sizeof(usint32_t), kmer_position_.nSizeIndex, fin),
              kmer_position_.nSizeIndex);
  FREAD_CHECK(fread(kmer_position_.index_proPos, sizeof(usint32_t), kmer_position_.nSizeIndex, fin),
              kmer_position_.nSizeIndex);

  fclose(fin);
  //TIME_INFO(TestIndex(), "Write Test Index");
}
