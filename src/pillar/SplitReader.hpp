/*
 * SplitReader.hpp
 *
 *  Created on: 2019. 10. 22.
 *      Author: Eun-Cheon Lim @ Postech Plant Genomics Lab.
 */

#ifndef PILLAR_SPLITREADER_HPP_
#define PILLAR_SPLITREADER_HPP_
#include <vector>
#include <string>
#include <boost/filesystem.hpp>
#include "../castle/TimeChecker.hpp"
#include "../castle/IOUtils.hpp"
#include "Kmer.hpp"
#include "SmallKmer.hpp"

namespace pillar {
using namespace std;
using namespace castle;
namespace bfs = boost::filesystem;
class SplitReader {
private:
	vector<string> paths;
	vector<uint64_t> file_sizes;
	vector<Type_Of_Format> file_types;
	vector<uint64_t> n_blocks;
	vector< vector<uint64_t> > skip_points;
	vector< vector<uint64_t> > actual_chunk_sizes;

	string working_path;
	uint32_t n_cores;
	int64_t n_dividing_elements;
	uint64_t total_file_size;
	uint64_t window_size;
	uint64_t minimizer_size;
public:
	SplitReader();
	~SplitReader();
	void split_a_read(const string& a_read, const string& a_qual);
	void split_a_read_alt(const string& a_read, const string& a_qual);
	SmallKmer get_single_minimizer(vector<SmallKmer>& kmers, const uint64_t start_id);
	SmallKmer get_minimizer(vector<SmallKmer>& kmers, const uint64_t start_id);
	Kmer get_minimizer_alt(vector<Kmer>& kmers, const uint64_t start_id);
	void add_read_path(const string& a_path);
	void set_working_path(const string& a_path);
	void prepare_parallel_processing();
	void calculate_all_file_sizes();
	void create_all_blocks();
	void create_splits();
	void find_a_split(const int64_t the_seq_id);
	void find_all_minimizers();
	void find_binning_minimizers(set<SmallKmer>& result_minimizers, const string& a_read);

	void print_progress(boost::mutex& a_mutex, int64_t& n_processing);
	template<class Ret>
	vector<Ret> process_a_block_task(const uint64_t the_seq_id, const string& func_name, const function<void(Ret&, boost::mutex&, string&, string&)>& a_callback);
};

template<class Ret>
vector<Ret> SplitReader::process_a_block_task(const uint64_t the_seq_id, const string& func_name, const function<void(Ret&, boost::mutex&, string&, string&)>& a_callback) {
	TimeChecker checker;
	string func_full_name = (boost::format("SplitReader.%s") % func_name).str();
	checker.setTarget(func_full_name);
	checker.start_without_output();
	vector<Ret> parallel_out;
	vector<function<void()> > tasks;
	string a_seq_path = paths[the_seq_id];
	const uint64_t max_block_index = actual_chunk_sizes[the_seq_id].size();
	cout << (boost::format("[%s] # blocks: %d\n") % func_full_name % max_block_index).str();
	cout << (boost::format("[%s] # seq: %s\n") % func_full_name % a_seq_path).str();
	boost::mutex print_mutex;
	boost::filesystem::path a_seq_local_path(a_seq_path);
	vector<uint64_t> out_n_succeeded_count(max_block_index + 1);
	vector<uint64_t> out_n_processed_count(max_block_index + 1);

	parallel_out.resize(max_block_index);
	set<uint64_t> debug_block;
	debug_block.insert(0);
	int64_t n_processing = 0;
	for (uint64_t block_id = 0; block_id < max_block_index; ++block_id) {
		tasks.push_back([&, &n_processing, the_seq_id, block_id] {
			bool debug = false;
			if(debug_block.end() == debug_block.find(block_id)) {
//				print_progress(print_mutex, n_processing);
				return;
			} else {
				debug = true;
			}
			auto& cur_out = parallel_out[block_id];
			const uint64_t skip_bytes = skip_points[the_seq_id][block_id];
			const uint64_t chunk_bytes = actual_chunk_sizes[the_seq_id][block_id];
			const auto a_file_type = file_types[the_seq_id];

			ifstream input_read_stream(a_seq_path, ios::binary);
			input_read_stream.seekg(skip_bytes, ios::beg);
//			vector<map<string, vector<int64_t>>> out_read;
			bfs::path a_local_path(a_seq_path);
//			string output_read_path = (boost::format("%s/%s.out.%d") % working_path % a_seq_local_path.filename().string() % block_id).str();
//			if(bfs::exists(output_read_path)) {
//				return;
//			}
//			ofstream out_read(output_read_path, ios::binary);
			string line;
			string read_id;

			uint64_t a_local_total_processed_bytes = 0;
			uint64_t a_local_processed_bytes = 0;
			stringstream s_read;
			stringstream s_qual;
			stringstream s_read_concatenating;
			stringstream s_qual_concatenating;

			uint64_t current_sequence_size = 0;
			uint64_t current_quality_size = 0;
			uint64_t state = 0;
			uint64_t local_processed_bytes = 0;

			uint64_t n_succeeded = 0;
			uint64_t n_processed_lines = 0;
			uint64_t n_processed_reads = 0;
			if(FASTQ == a_file_type) {
				while(getline(input_read_stream, line, '\n')) {
					local_processed_bytes += line.size() + 1;
					++n_processed_lines;
					a_local_total_processed_bytes += line.size() + 1;
					if ('@' == line[0] && current_sequence_size == current_quality_size) {
						state = 1;
						if (0 == current_sequence_size) {
							continue;
						}
						++n_processed_reads;
//						if(debug) {
//							if (0 == (n_processed_reads & 1023)) {
//								cout << (boost::format("[%s] # processed: %d\n") % func_full_name % n_processed_reads).str();
//							}
//						}
						string a_read = s_read_concatenating.str();
						IOUtils::replace_unknown_bases(a_read);
						string a_qual = s_qual_concatenating.str();
						int64_t base_start_id = 0;
						int64_t base_end_id = a_read.size() - 1;
						int64_t new_read_length = IOUtils::get_trimmed_indexes(base_start_id, base_end_id, a_read);
						s_read.write(&a_read[base_start_id], new_read_length);
						s_qual.write(&a_qual[base_start_id], new_read_length);
						a_local_processed_bytes += (new_read_length << 1) + 2;
						s_read << "\n";
						s_qual << "\n";


						if (a_local_processed_bytes > IOUtils::PARTIAL_READ_BUFFER_SIZE) {
							s_read.str(string());
							s_qual.str(string());
							a_local_processed_bytes = 0;
						}
						a_callback(cur_out, print_mutex, a_read, a_qual);

						s_read_concatenating.str(string());
						s_qual_concatenating.str(string());
						current_sequence_size = 0;
						current_quality_size = 0;
					} else if ('+' == line[0] && 1 == state) {
						state = 2;
					} else if (1 == state) {
						current_sequence_size += line.size();
						s_read_concatenating << line;
					} else if (2 == state) {
						current_quality_size += line.size();
						s_qual_concatenating << line;
					}
					if(a_local_total_processed_bytes >= chunk_bytes) {
						break;
					}
				}
				string a_read = s_read_concatenating.str();
				if (a_read.size() > 0) {
					IOUtils::replace_unknown_bases(a_read);
					string a_qual = s_qual_concatenating.str();
					++n_processed_reads;
					a_callback(cur_out, print_mutex, a_read, a_qual);
				}
			} else if (FASTA == a_file_type) {
				while(getline(input_read_stream, line, '\n')) {
					local_processed_bytes += line.size() + 1;
					a_local_total_processed_bytes += line.size() + 1;
					++n_processed_lines;
					if ('>' == line[0]) {
						++n_processed_reads;
//						if(debug) {
//							if (0 == (n_processed_reads & 1023)) {
//								cout << (boost::format("[%s] # processed: %d\n") % func_full_name % n_processed_reads).str();
//							}
//						}
						if(0 == current_sequence_size) {
							read_id = line;
							continue;
						}
						string a_read = s_read_concatenating.str();
						string a_qual;
						a_callback(cur_out, print_mutex, a_read, a_qual);

						s_read_concatenating.str(string());
						current_sequence_size = 0;
						read_id = line;
					} else {
						current_sequence_size += line.size();
						s_read_concatenating << line;
					}

					if(a_local_total_processed_bytes >= chunk_bytes) {
						break;
					}
				}
				string a_read = s_read_concatenating.str();
				if(a_read.size() > 0) {
					++n_processed_reads;
					string a_qual;
					a_callback(cur_out, print_mutex, a_read, a_qual);
				}
			}

			out_n_succeeded_count[block_id] = n_succeeded;
			out_n_processed_count[block_id] = n_processed_reads;
			print_progress(print_mutex, n_processing);
		});
	}
	ParallelRunner::run_unbalanced_load(n_cores, tasks);
	return parallel_out;
}

} /* namespace pillar */

#endif /* PILLAR_SPLITREADER_HPP_ */
