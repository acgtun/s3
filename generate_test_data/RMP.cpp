/*
 * RMP (Reads Mapping Protein Database)
 * Haifeng Chen
 * University of Southern California
 * August 4, 2014
 *
 */

/*
 * File name: all the letters are lower case
 * Function name: the first letter of every word are upper case, \
 *   the others are lower case
 * Variable name: the first letter of each word except \
 *   the first word are upper case, the others are lower case
 *
 * make it as simple as possible
 */

#include "option.h"
#include "buildindex.h"
#include "dbsearch.h"

int main(int argc, const char* argv[]) {
	/* (1) Get parameters */
	Option opt(argc, argv);
	if(opt.opt_error) return 0;

	Index index(opt);

	if (opt.bIndexExist) {
		/* (2) Read Reference and Hash Table from index */
		TIME_INFO(index.ReadIndex(), "Read Index");
	} else {
		/* (3) Read Reference */
		TIME_INFO(index.BuildIndex(), "Read Protein Database and Build Index");
	}
	if(argc <= 2) return 0;
	/* (5) Build Reference Index */
	DBSearch dbsearch(opt, &index);
	TIME_INFO(dbsearch.Search(), "Search");

	/* release memory*/
	//free(proDB.proteinSeq);
	//free(hashTable.counter);
	//free(hashTable.index);
	return 0;
}
