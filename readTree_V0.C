// V0 reconstruction macro from local aurora trees
// OliverM 2018 Lund

#include <iostream>

using namespace std;

void readTree_V0(Int_t nEvents=10, const Char_t *inputFile="test.list", const Char_t *outputFile="test.root") {

	gROOT->LoadMacro("$HOME/sq/load_libraries.C");
	load_libraries();

	printf(" WHAT IS UP \n", );
}