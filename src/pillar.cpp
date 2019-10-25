//============================================================================
// Name        : pillar.cpp
// Author      : Eun-Cheon Lim @ Postech Plant Genomics Lab.
// Version     :
// Copyright   : All rights reserved to Eun-Cheon Lim
// Description : Hello World in C++, Ansi-style
//============================================================================

#include "pillar/SplitReader.hpp"
#include <iostream>
using namespace pillar;
using namespace std;

int main(int argc, char **argv) {
	setbuf(stdout, NULL);
	TimeChecker checker;
	checker.setTarget("Chronos.Main");
	checker.start();
	string input_path("/home/pgr/data_1/LEC/genomics/asmb/Ler-0/Ler-0.fasta");
	string a_working_path("/data_1/LEC/eclipse_cdt/workspace/pillar/src/test");
	SplitReader sr;
	sr.add_read_path(input_path);
	sr.set_working_path(a_working_path);
	sr.prepare_parallel_processing();
	sr.find_all_minimizers();
////	sr.create_splits();
	cout << checker;
	return 0;
}
