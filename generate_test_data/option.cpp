#include "option.h"

void Option::PrintSynopsis() {
	printf(
			"The input command is incorrect.\nFor more info, please check: http://code.google.com/p/perm/\n");
}

int Option::GetIntVal(int argc, const char** argv, const char * str,
		int & nVal) {
	for (int i = 1; i < argc; i++) {
		if (strcmp(argv[i], str) == 0 && i + 1 <= argc - 1) {
			if (argv[i + 1][0] != '-') {
				nVal = atoi(argv[i + 1]);
				return 1;
			}
		}
	}
	return 0;
}

int Option::GetStrVal(int argc, const char** argv, const char* str,
		string & strVal) {
	for (int i = 1; i < argc; i++) {
		if (strcmp(argv[i], str) == 0 && i + 1 <= argc - 1) {
			if (argv[i + 1][0] != '-') {
				strVal = argv[i + 1];
				return 1;
			}
		}
	}
	return 0;
}

int Option::ChkStrExist(int argc, const char** argv, const char* str) {
	for (int i = 1; i < argc; i++) {
		if (strcmp(argv[i], str) == 0) {
			return 1;
		}
	}
	return 0;
}

void Option::GetNameBeforeDot(const string & strFile, string & fileName) {
	int pos = strFile.find_last_of('.');
	fileName = strFile.substr(0, pos);
}

void Option::GetNumOfReads() {
	ifstream fin(readsFile.c_str());
	if (!fin.good()) {
		printf("--ERROR INFO-- reads file open error. %s\n", readsFile.c_str());
		exit(EXIT_FAILURE);
	}
	nNumOfreads = 0;
	string strLine;
	while (getline(fin, strLine)) {
		if (strLine[0] == '>')
			nNumOfreads++;
	}
	fin.close();
	INFO("There are", nNumOfreads, "reads in the reads file");
}

void Option::PrintOutVersion() {
	printf("####################################################\n");
	printf("#                                                  #\n");
	printf("#  RMP (Reads Mapping Protein Database)            #\n");
	printf("#  Haifeng Chen (haifengc at usc dot edu)          #\n");
	printf("#  University of Southern California               #\n");
	printf("#  Aug 24, 2014                                    #\n");
	printf("#                                                  #\n");
	printf("####################################################\n");
	printf("\n");
}

Option::Option(int argc, const char* argv[]) {
	opt_error = false;
	if (argc <= 1 || argv == NULL) {
		PrintSynopsis();
		opt_error = true;
	} else {
		PrintOutVersion();
		cout << "--INFO-- Input command:";
		for (int i = 0; i < argc; i++) {
			cout << " " << argv[i];
		}
		cout << endl;

		refFile = argv[1];
		cout << "--INFO-- The reference file is " << refFile << endl;

		bIndexExist = 0;
		string indexsf = ".sdindex";
		if (refFile.size() > indexsf.size()
				&& refFile.substr(refFile.size() - indexsf.size()) == indexsf) {
			bIndexExist = 1;
		}

		bSaveIndex = 0;
		if (ChkStrExist(argc, argv, "-s")) {
			bSaveIndex = 1;
		}

		bSaveIndex = 1;

		if (GetStrVal(argc, argv, "-s", indexFile) == 0 && bSaveIndex == 1) {
			string fileName;
			GetNameBeforeDot(refFile, fileName);
			indexFile = fileName;
			indexFile += indexsf;
			//cout << "indexFile = " << indexFile << endl;
		}

		if (GetIntVal(argc, argv, "-e", expected_hits) == 0) {
			expected_hits = 50;
			//printf("TEST The expected number of his for each query is %d.\n", expected_hits);
			//ERROR_INFO("please input the number of expected hits for each query");
		}
		INFO("The expected number of hits for each query is", expected_hits);

		if (GetStrVal(argc, argv, "-q", readsFile)) {
			cout << "--INFO-- The reads file is " << readsFile << "." << endl;

			if (GetStrVal(argc, argv, "-o", outputFile) == 0) {
				string fileName;
				GetNameBeforeDot(readsFile, fileName);
				outputFile = fileName;
				outputFile += ".sam";
			}
			GetNumOfReads();
		}
	}
}

Option::~Option() {
}
