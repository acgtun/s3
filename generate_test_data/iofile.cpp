#include "iofile.h"

usint64_t GetLineFromString(const char * strVal, char * strRet) {
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

usint64_t ReadWholeFile(const string & fileName, char **strRet) {
	FILE * fin = fopen(fileName.c_str(), "rb");
	FILE_OPEN_CHECK(fin);
	fseek(fin, 0, SEEK_END);
	usint64_t size = ftell(fin);
	MEMORY_ALLOCATE_CHECK(*strRet = (char* ) malloc(sizeof(char) * (size + 1)));
	fseek(fin, 0, SEEK_SET);
	FREAD_CHECK(fread(*strRet, 1, size, fin), size);
	fclose(fin);
	return size;
}

