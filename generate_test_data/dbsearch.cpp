#include "dbsearch.h"
#include "iofile.h"
#include <stdlib.h>

DBSearch::DBSearch(const Option& opt, const Index * index) :
		opt_(opt), index_(index) {
	reads_.nReadsNum = opt_.nNumOfreads;
	ReadReads();
	for (usint32_t i = 0; i < MAX_QUERY_LEN / 6 + 1; i++) {
		seed_start_pos_[i] = i * HASHAALEN;
	}
}

DBSearch::~DBSearch() {
	for (usint32_t i = 0; i < reads_.nReadsNum; i++) {
		free(reads_.readsArray[i]);
	}
	free(reads_.readsArray);
}

void DBSearch::ReadReads() {
	char * strReads;
	usint64_t readsLen = ReadWholeFile(opt_.readsFile, &strReads);
	MEMORY_ALLOCATE_CHECK(reads_.readsArray = (char ** ) malloc(reads_.nReadsNum * sizeof(char *)));

	INFO("Read Reads from", opt_.readsFile);
	usint32_t readID = 0, readLen;
	char strRead[MAX_LINE_LEN];
	string strReadApp;
	for (usint64_t i = 0; i < readsLen; i++) {
		readLen = GetLineFromString(&strReads[i], strRead);
		i += readLen;
		if (readLen == 0)
			continue;
		if (strRead[0] == '>') {
			if (strReadApp.size() != 0) {
				MEMORY_ALLOCATE_CHECK(
						reads_.readsArray[readID] = (char * ) malloc((strReadApp.size() + 1) * sizeof(char)));
				memcpy(reads_.readsArray[readID], strRead, strReadApp.size());
				reads_.readsArray[readID][strReadApp.size()] = 0;
				reads_.vReadSize[readID] = strReadApp.size();
				readID++;
				strReadApp.clear();
			}
			continue;
		} else {
			strReadApp += strRead;
		}
	}
	if (strReadApp.size() != 0) {
		MEMORY_ALLOCATE_CHECK(reads_.readsArray[readID] = (char * ) malloc((strReadApp.size() + 1) * sizeof(char)));
		memcpy(reads_.readsArray[readID], strRead, strReadApp.size());
		reads_.readsArray[readID][strReadApp.size()] = 0;
		reads_.vReadSize[readID] = strReadApp.size();
		readID++;
	}
	free(strReads);
}

bool SortCMP(const pair<usint32_t, usint32_t> & a, const pair<usint32_t, usint32_t> & b) {
	if (a.first == b.first) {
		return a.second < b.second;
	}
	return a.first < b.first;
}

void DBSearch::Search() {
	ofstream fout("results.txt");
	for (usint32_t i = 0; i < reads_.nReadsNum; i++) {
		usint32_t hashValue = 0, numSeed = reads_.vReadSize[i] / HASHAALEN;
		vector<pair<usint32_t, usint32_t> > results;
		fout << "read=" << reads_.readsArray[i] << endl;
		for (usint32_t j = 0; j < numSeed; j++) {
			hashValue = GetHashValue(&(reads_.readsArray[i][seed_start_pos_[j]]));
			for (usint64_t c = index_->kmer_nearest_.counter[hashValue];
					c < index_->kmer_nearest_.counter[hashValue + 1]; c++) {
				usint32_t Kmer = index_->kmer_nearest_.index[c];

				usint32_t nearestKmerS = index_->kmer_position_.counter[Kmer];
				usint32_t nearestKmerE = index_->kmer_position_.counter[Kmer + 1];
				for (usint64_t c2 = nearestKmerS; c2 < nearestKmerE; c2++) {
					results.push_back(
							make_pair(index_->kmer_position_.index_proID[c2], index_->kmer_position_.index_proPos[c2]));
				}
			}
		}
		sort(results.begin(), results.end(), SortCMP);

		if (results.size() > 0) {
			fout << "protein ID: " << results[0].first << endl;
			fout << results[0].second;
		}
		for (usint32_t i = 1; i < results.size(); i++) {
			if (results[i].first == results[i - 1].first) {
				fout << " " << results[i].second;
			} else {
				fout << endl;
				fout << "protein ID: " << results[i].first << endl;
				fout << results[i].second;
			}
		}
		fout << endl;
		for (vector<pair<usint32_t, usint32_t> >::iterator it = results.begin(); it != results.end(); it++) {
			fout << it->first << " " << it->second << " ";
			for (usint32_t j = it->second; j < index_->proDB_.vProSize[it->first]; j++) {
				fout << index_->proDB_.vProSeq[it->first][j];
			}
			fout << " " << index_->proDB_.vProName[it->first] << endl;
		}
	}
	fout.close();
}

