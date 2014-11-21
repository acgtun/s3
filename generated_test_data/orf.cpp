#include "orf.h"
#include <map>

void SetGeneticCodes(map<string, char> & mapGeneticCodes) {
	char code[4];
	for (usint32_t i = 0; i < Base1.size(); i++) {
		code[0] = Base1[i];
		code[1] = Base2[i];
		code[2] = Base3[i];
		code[3] = 0;
		mapGeneticCodes[code] = AAs[i];
	}
}

char complimentBase(char nt) {
	switch (nt) {
	case 'A':
		return ('T');
	case 'C':
		return ('G');
	case 'G':
		return ('C');
	case 'T':
		return ('A');
	default:
		ERROR_INFO("DNA has non-ACGT characters");
		return 'N';
	}
}

void ReverseCompliment(string & strRead_rev, const string & strRead, const int & len) {
	for (int i = 0; i < len; i++) {
		strRead_rev += complimentBase(strRead[len - i - 1]);
	}
}

void orf6(const string & querySeq, map<string, char> & mapGeneticCodes, vector<string> & translatedAA) {
	cout << "querySeq = " << querySeq << endl;
	int len = querySeq.size();
	string strAA;
	char chrAA;
	len -= 3;
	for (int s = 0; s < 3; s++) {
		strAA.clear();
		for (int i = s; i <= len; i += 3) {
			chrAA = mapGeneticCodes[querySeq.substr(i, 3)];
			if (chrAA == '*')
				break;
			strAA += chrAA;
		}
		translatedAA[s] = strAA;
	}

	string querySeq_rev;
	ReverseCompliment(querySeq_rev, querySeq, querySeq.size());

	for (int s = 0; s < 3; s++) {
		strAA.clear();
		for (int i = s; i <= len; i += 3) {
			chrAA = mapGeneticCodes[querySeq_rev.substr(i, 3)];
			if (chrAA == '*')
				break;
			strAA += chrAA;
		}
		translatedAA[s + 3] = strAA;
	}
}
