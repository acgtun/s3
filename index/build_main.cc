#include "./../util/option.h"
#include "index.h"

int main(int argc, const char* argv[]) {
  InitProgram(argc, argv);

  Index index;

  TIME_INFO(index.BuildIndex(), "Read Protein Database and Build Index");

  return 0;
}
