/*
 * MinimerSorter.hpp
 *
 *  Created on: Sep 24, 2019
 *      Author: Eun-Cheon Lim @ Postech Plant Genomic Recombination Lab.
 *      Email : abysslover@gmail.com
 *    License : GPLv3
 */

#ifndef MINIMERSORTER_HPP_
#define MINIMERSORTER_HPP_
#include <string>
#include <vector>
#include <stack>
#include <queue>
#include <iostream>
#include <functional>
#include <memory>
#include <bitset>
#include <parallel/algorithm>
#include <boost/filesystem.hpp>
#include <boost/thread/thread.hpp>
#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>
#include <read_sorter/castle/BinaryDecoder.hpp>
#include <read_sorter/castle/IOUtils.hpp>
#include <read_sorter/castle/StringUtils.hpp>
#include <read_sorter/castle/TimeChecker.hpp>

namespace pillar {
	using namespace castle;
	namespace bfs = boost::filesystem;
	using namespace std;

class MinimerSorter {
private:
		uint64_t k;
		uint32_t n_cores;
		uint64_t empty_kmer;
		int64_t n_dividing_elements;
		uint64_t total_file_size;

		string working_path;
		vector<uint64_t> file_sizes;
		vector<Type_Of_Format> file_types;
		vector<string> paths;
		vector<uint64_t> n_blocks;
		vector< vector<uint64_t> > skip_points;
		vector< vector<uint64_t> > actual_chunk_sizes;
public:
	MinimerSorter();
	~MinimerSorter();
	void set_k(const uint64_t a_k);
	void set_working_path(const string& a_path);
	void add_read_path(const string& a_path);
	void prepare_parallel_processing();
	void find_premer_connections();
	void connect_premers(PremerMap& a_premer_connection_map, const string& a_read, BinaryEncoder& be);
	void find_a_premer_connection(const int64_t the_seq_id);
	void sort_all_reads();
	void sort_a_read(const int64_t the_seq_id);
	void merge_sorted_files(const int64_t the_seq_id);
	bool process_a_fastq_line(ofstream& out, boost::mutex& a_mutex, const string& a_read, const string& a_qual);
	bool process_a_fasta_line(ofstream& out, boost::mutex& a_mutex, const string& a_read);
	void calculate_all_file_sizes();
	void create_all_blocks();
	void print_progress(boost::mutex& a_mutex, int64_t& n_processing);
};

} /* namespace pillar */

#endif /* MINIMERSORTER_HPP_ */
