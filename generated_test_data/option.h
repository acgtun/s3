#pragma once
#ifndef PARAMETER_H_
#define PARAMETER_H_

#include "sdk.h"

class Option {
public:
	bool opt_error;
	//input and out parameter
	string readsFile;
	string refFile;

	usint64_t nNumOfreads;

	string outputFile;
	string indexFile;
	int bSaveIndex; // Default is false
	int bIndexExist;

	int expected_hits;

	Option(int argc, const char* argv[]);
	~Option();

private:
	void PrintSynopsis();
	int GetIntVal(int argc, const char** argv, const char * str, int & nVal);
	int GetStrVal(int argc, const char** argv, const char* str, string & strVal);
	int ChkStrExist(int argc, const char** argv, const char* str);
	void GetNameBeforeDot(const string & strFile, string & fileName);
	void GetNumOfReads();
	void PrintOutVersion();
	void GetThresholdScore();
};

#endif /* PARAMETER_H_ */
