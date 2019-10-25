//============================================================================
// Name        : pillar.cpp
// Author      : Eun-Cheon Lim @ Postech Plant Genomics Lab.
// Version     :
// Copyright   : Copyright 2019- All rights reserved.
// Description : Hello World in C++, Ansi-style
//============================================================================

#include "pillar/MinimerSorter.hpp"

using namespace std;
using namespace pillar;
using namespace castle;
int main(int argc, char **argv) {
	setbuf(stdout, NULL);
	TimeChecker checker;
	checker.setTarget("Chronos.Main");
	checker.start();
	string input_path("/home/vincent/data_1/genomics/asmb/Ler-0/Ler-0.fasta");
	string a_working_path("/home/vincent/data_1/genomics/asmb/Ler-0");
	MinimerSorter ms;
	ms.set_k(15);
	ms.set_working_path(a_working_path);
	ms.add_read_path(input_path);
	ms.prepare_parallel_processing();
	ms.find_premer_connections();
//	ms.sort_all_reads();
	cout << checker;
}
