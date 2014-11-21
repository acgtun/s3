#pragma once
#ifndef IOFILE_H_
#define IOFILE_H_

#include "sdk.h"

usint64_t ReadWholeFile(const string & fileName, char **strRet);
usint64_t GetLineFromString(const char * strVal, char * strRet);

#endif /* IOFILE_H_ */
