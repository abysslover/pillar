/*
 * MinimerSorter.cpp
 *
 *  Created on: Sep 24, 2019
 *      Author: Eun-Cheon Lim @ Postech Plant Genomic Recombination Lab.
 *      Email : abysslover@gmail.com
 *    License : GPLv3
 */

#include "MinimerSorter.hpp"

namespace pillar {

MinimerSorter::MinimerSorter() : k(15), n_dividing_elements(100), total_file_size(0) {
	TimeChecker checker;
	n_cores = checker.get_number_of_cores();
	BinaryEncoder kmer_encoder;
	kmer_encoder.setReadLength(32);
	empty_kmer = kmer_encoder.encodeToTwoBitsLong("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAC", 0);
}

MinimerSorter::~MinimerSorter() {
}

void MinimerSorter::set_k(const uint64_t a_k) {
	k = a_k;
}
void MinimerSorter::set_working_path(const string& a_path) {
	working_path = a_path;
	if (!bfs::exists(working_path)) {
		bfs::create_directories(working_path);
	}
}
void MinimerSorter::add_read_path(const string& a_path) {
	if (!bfs::exists(a_path)) {
		return;
	}
	paths.push_back(a_path);
}

void MinimerSorter::prepare_parallel_processing() {
	calculate_all_file_sizes();
	create_all_blocks();
}
void MinimerSorter::find_premer_connections() {
	for(uint64_t a_seq_id = 0; a_seq_id < paths.size(); ++ a_seq_id) {
		find_a_premer_connection(a_seq_id);
	}
}
void MinimerSorter::connect_premers(PremerMap& a_premer_connection_map, const string& a_read, BinaryEncoder& be) {
	const int64_t a_step = 400;
	int64_t base_start_pos = 0;
	int64_t base_end_pos = a_step;
	int64_t actual_base_end_pos = a_step;
	int64_t base_max_pos = a_read.size();
	vector<uint64_t> serial_connections_forward;
	vector<uint64_t> serial_connections_reverse;
	vector<uint64_t> results;
//	cout << (boost::format("[MinimerSorter.connect_premers] Length of a query: %d\n") % a_read.size());
	while(base_end_pos < base_max_pos) {
		base_end_pos = base_start_pos + a_step;
		actual_base_end_pos = min(base_end_pos, base_max_pos);
//		cout << (boost::format("[MinimerSorter.connect_premers] Range of the query: %d~%d\n") % base_start_pos % actual_base_end_pos);
		be.get_all_for_rev_premers(a_read, base_start_pos, actual_base_end_pos, results);
		serial_connections_forward.push_back(results[0]);
		serial_connections_reverse.push_back(results[1]);
		base_start_pos = base_end_pos;
	}
	int64_t max_pillar_id = serial_connections_forward.size();

//	BinaryDecoder bd;
//	bd.set_k(be.getReadLength());
	for(int64_t pillar_id = 0; pillar_id < max_pillar_id - 1; ++pillar_id) {
		auto& cur_premer = serial_connections_forward[pillar_id];
		auto& next_premer = serial_connections_forward[pillar_id + 1];
		auto& cur_map_container = a_premer_connection_map[cur_premer];
//		if(cur_map_container.smallest_key > cur_premer) {
//			cur_map_container.smallest_key = cur_premer;
//		}
//		if(cur_map_container.largest_key < cur_premer) {
//			cur_map_container.largest_key = cur_premer;
//		}
		++cur_map_container.next[next_premer];
		auto& next_map_container = a_premer_connection_map[next_premer];
		++next_map_container.prev[cur_premer];
	}
	for(int64_t pillar_id = max_pillar_id - 1; pillar_id > 0; --pillar_id) {
		auto& cur_premer = serial_connections_reverse[pillar_id];
		auto& next_premer = serial_connections_reverse[pillar_id - 1];
		auto& cur_map_container = a_premer_connection_map[cur_premer];
//		if(cur_map_container.smallest_key > cur_premer) {
//			cur_map_container.smallest_key = cur_premer;
//		}
//		if(cur_map_container.largest_key < cur_premer) {
//			cur_map_container.largest_key = cur_premer;
//		}
		++cur_map_container.next[next_premer];
		auto& next_map_container = a_premer_connection_map[next_premer];
		++next_map_container.prev[cur_premer];
	}
}
void MinimerSorter::find_a_premer_connection(const int64_t the_seq_id) {
	TimeChecker checker;
	checker.setTarget("MinimerSorter.find_a_premer_connection");
	checker.start_without_output();
	vector<function<void()> > tasks;
	string a_seq_path = paths[the_seq_id];
	const uint64_t max_block_index = actual_chunk_sizes[the_seq_id].size();
	cout << (boost::format("[MinimerSorter.find_a_premer_connection] # blocks: %d\n") % max_block_index).str();
	boost::mutex print_mutex;
	boost::filesystem::path a_seq_local_path(a_seq_path);
	vector<uint64_t> out_n_succeeded_count(max_block_index + 1);
	vector<uint64_t> out_n_processed_count(max_block_index + 1);

//	set<uint64_t> debug_block;
//	debug_block.insert(0);
//	debug_block.insert(106);
//	debug_block.insert(319);
	BinaryEncoder be;
	be.setReadLength(k);
	BinaryDecoder bd;
	bd.set_k(k);
	int64_t n_processing = 0;
	for (uint64_t block_id = 0; block_id < max_block_index; ++block_id) {
		tasks.push_back([&, &n_processing, the_seq_id, block_id] {
//			if(debug_block.end() == debug_block.find(block_id)) {
//				return;
//			}
			const uint64_t skip_bytes = skip_points[the_seq_id][block_id];
			const uint64_t chunk_bytes = actual_chunk_sizes[the_seq_id][block_id];
			const auto a_file_type = file_types[the_seq_id];

			ifstream input_read_stream(a_seq_path, ios::binary);
			input_read_stream.seekg(skip_bytes, ios::beg);
//			vector<map<string, vector<int64_t>>> out_read;
			bfs::path a_local_path(a_seq_path);
			string output_read_path = (boost::format("%s/%s.out.%d") % working_path % a_seq_local_path.filename().string() % block_id).str();
//			if(bfs::exists(output_read_path)) {
//				return;
//			}
			ofstream out_read(output_read_path, ios::binary);
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
			PremerMap premer_connection_map;

			premer_connection_map.set_empty_key(empty_kmer);
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

						connect_premers(premer_connection_map, a_read, be);

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
					connect_premers(premer_connection_map, a_read, be);
				}
			} else if (FASTA == a_file_type) {
				while(getline(input_read_stream, line, '\n')) {
					local_processed_bytes += line.size() + 1;
					a_local_total_processed_bytes += line.size() + 1;
					++n_processed_lines;
					if ('>' == line[0]) {
						++n_processed_reads;
						if(0 == current_sequence_size) {
							read_id = line;
							continue;
						}
						string a_read = s_read_concatenating.str();
						connect_premers(premer_connection_map, a_read, be);

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
					connect_premers(premer_connection_map, a_read, be);
				}
			}

			out_n_succeeded_count[block_id] = n_succeeded;
			out_n_processed_count[block_id] = n_processed_reads;
//			{
//				boost::lock_guard<boost::mutex> lock(print_mutex);
//				cout << (boost::format("[BlockReader.process_a_fastq_file_to_file] DONE: %d-%d: %d~%d(%d)\n")
//				% the_seq_id % block_id % skip_bytes % (skip_bytes + chunk_bytes) % chunk_bytes).str();
//			}
//			sort(chunk_premers.begin(), chunk_premers.end());
//			__gnu_parallel::sort(chunk_premers.begin(), chunk_premers.end());
			for(auto& a_premer_pair: premer_connection_map) {
				stringstream ss;

				stringstream ss_prev;
				for(auto& a_next_premer_pair: a_premer_pair.second.prev) {
					ss_prev << (boost::format(":%s_%d") % bd.decode(a_next_premer_pair.first) % a_next_premer_pair.second).str();
				}
				string prev_str = ss_prev.str();
				if(prev_str.size() > 0) {
					ss << prev_str.substr(1);
				}
				ss << "\t" << bd.decode(a_premer_pair.first) << "\t";

				stringstream ss_next;

				for(auto& a_next_premer_pair: a_premer_pair.second.next) {
					ss_next << (boost::format(":%s_%d") % bd.decode(a_next_premer_pair.first) % a_next_premer_pair.second).str();
				}
				string next_str = ss_next.str();
				if(next_str.size() > 0) {
					ss << next_str.substr(1);
				}
				ss << "\n";
				out_read << ss.str();
			}
			print_progress(print_mutex, n_processing);
		});
	}
	ParallelRunner::run_unbalanced_load(n_cores, tasks);
//	ParallelRunner::run_unbalanced_load(1, tasks);
//	ParallelRunner::run_step_wise(1, tasks, ParallelRunner::empty_ranged_func);

//	uint64_t n_total_succeeded = 0;
//	uint64_t n_total_processed = 0;
//	for (uint64_t block_id = 0; block_id < max_block_index; ++block_id) {
//		uint64_t n_cur_succeeded = out_n_succeeded_count[block_id];
//		uint64_t n_cur_processed = out_n_processed_count[block_id];
//		n_total_succeeded += n_cur_succeeded;
//		n_total_processed += n_cur_processed;
//	}
//	cout << (boost::format("[BlockReader.find_a_premer_connection] # reads: %d\n") % n_total_processed).str();
//	string output_suffix = (boost::format("%s/%s.out") % working_path % a_seq_local_path.filename().string()).str();
//	merge_sorted_files(the_seq_id);
}
void MinimerSorter::sort_all_reads() {
	for(uint64_t a_seq_id = 0; a_seq_id < paths.size(); ++ a_seq_id) {
		sort_a_read(a_seq_id);
	}
}

// use this code as a template
void MinimerSorter::sort_a_read(const int64_t the_seq_id) {
	TimeChecker checker;
	checker.setTarget("MinimerSorter.sort_a_read");
	checker.start_without_output();
	vector<function<void()> > tasks;
	string a_seq_path = paths[the_seq_id];
	const uint64_t max_block_index = actual_chunk_sizes[the_seq_id].size();
	cout << (boost::format("[MinimerSorter.sort_a_read] # blocks: %d\n") % max_block_index).str();
	boost::mutex print_mutex;
	boost::filesystem::path a_seq_local_path(a_seq_path);
	vector<uint64_t> out_n_succeeded_count(max_block_index + 1);
	vector<uint64_t> out_n_processed_count(max_block_index + 1);

//	set<uint64_t> debug_block;
//	debug_block.insert(0);
//	debug_block.insert(106);
//	debug_block.insert(319);
	BinaryEncoder be;
	be.setReadLength(k);

	for (uint64_t block_id = 0; block_id < max_block_index; ++block_id) {
		tasks.push_back([&, the_seq_id, block_id] {
//			if(debug_block.end() == debug_block.find(block_id)) {
//				return;
//			}
			const uint64_t skip_bytes = skip_points[the_seq_id][block_id];
			const uint64_t chunk_bytes = actual_chunk_sizes[the_seq_id][block_id];
			const auto a_file_type = file_types[the_seq_id];

			ifstream input_read_stream(a_seq_path, ios::binary);
			input_read_stream.seekg(skip_bytes, ios::beg);
//			vector<map<string, vector<int64_t>>> out_read;
			bfs::path a_local_path(a_seq_path);
			string output_read_path = (boost::format("%s/%s.out.%d") % working_path % a_seq_local_path.filename().string() % block_id).str();
//			if(bfs::exists(output_read_path)) {
//				return;
//			}
			ofstream out_read(output_read_path, ios::binary);
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
			vector<Premer<string>> chunk_premers;
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

						chunk_premers.push_back(be.get_premer(a_read, 0, a_read.size(), 'C'));
//						if(process_a_fastq_line(out_read, print_mutex, a_read, a_qual)) {
//							++n_succeeded;
//						}
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
					chunk_premers.push_back(be.get_premer(a_read, 0, a_read.size(), 'C'));
				}
			} else if (FASTA == a_file_type) {
				while(getline(input_read_stream, line, '\n')) {
					local_processed_bytes += line.size() + 1;
					a_local_total_processed_bytes += line.size() + 1;
					++n_processed_lines;
					if ('>' == line[0]) {
						++n_processed_reads;
						if(0 == current_sequence_size) {
							read_id = line;
							continue;
						}
						string a_read = s_read_concatenating.str();
//						if(process_a_fasta_line(out_read, print_mutex, a_read)) {
//							++n_succeeded;
//						}
						chunk_premers.push_back(be.get_premer(a_read, 0, a_read.size(), 'C'));
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
					chunk_premers.push_back(be.get_premer(a_read, 0, a_read.size(), 'C'));
				}
			}

			out_n_succeeded_count[block_id] = n_succeeded;
			out_n_processed_count[block_id] = n_processed_reads;;
//			{
//				boost::lock_guard<boost::mutex> lock(print_mutex);
//				cout << (boost::format("[BlockReader.process_a_fastq_file_to_file] DONE: %d-%d: %d~%d(%d)\n")
//				% the_seq_id % block_id % skip_bytes % (skip_bytes + chunk_bytes) % chunk_bytes).str();
//			}
//			sort(chunk_premers.begin(), chunk_premers.end());
			__gnu_parallel::sort(chunk_premers.begin(), chunk_premers.end());
			for(auto& a_premer: chunk_premers) {
				out_read << a_premer.m_key << "\t" << a_premer.n_pos << "\t" << a_premer.m_val << "\n";
//				out_read << bd.decode(a_premer.m_key) << "\t" << a_premer.n_pos << "\t" << a_premer.m_val << "\n";
			}
		});
	}
	ParallelRunner::run_unbalanced_load(n_cores, tasks);
//	ParallelRunner::run_unbalanced_load(1, tasks);
//	ParallelRunner::run_step_wise(1, tasks, ParallelRunner::empty_ranged_func);

	uint64_t n_total_succeeded = 0;
	uint64_t n_total_processed = 0;
	for (uint64_t block_id = 0; block_id < max_block_index; ++block_id) {
		uint64_t n_cur_succeeded = out_n_succeeded_count[block_id];
		uint64_t n_cur_processed = out_n_processed_count[block_id];
		n_total_succeeded += n_cur_succeeded;
		n_total_processed += n_cur_processed;
	}
	cout << (boost::format("[MinimerSorter.sort_a_read] # reads: %d\n") % n_total_processed).str();
	string output_suffix = (boost::format("%s/%s.out") % working_path % a_seq_local_path.filename().string()).str();
	merge_sorted_files(the_seq_id);
//	IOUtils::combine_intermediate_files(output_suffix, max_block_index);
//	IOUtils::remove_intermediate_files(output_suffix, max_block_index, n_cores);
}

static bool ComparePremerPair(pair<Premer<string>, uint64_t> lhs, pair<Premer<string>, uint64_t> rhs) {
    return lhs.first < rhs.first;
}
void MinimerSorter::merge_sorted_files(const int64_t the_seq_id) {
	TimeChecker checker;
	checker.setTarget("BlockReader.merge_sorted_files");
	checker.start_without_output();
	string a_seq_path = paths[the_seq_id];
	bfs::path a_seq_local_path(a_seq_path);
	const uint64_t max_block_index = actual_chunk_sizes[the_seq_id].size();
    priority_queue<pair<Premer<string>, uint64_t>, vector<pair<Premer<string>, uint64_t>>, std::function<bool(pair<Premer<string>, uint64_t>, pair<Premer<string>, uint64_t>)>> min_heap(ComparePremerPair);
    ifstream* input_files = new ifstream[max_block_index];
    for (uint64_t block_id = 0; block_id < max_block_index; ++block_id) {
		string output_read_path = (boost::format("%s/%s.out.%d") % working_path % a_seq_local_path.filename().string() % block_id).str();
		input_files[block_id].open(output_read_path.c_str());
		Premer<string> min_premer_in_this_block;
		input_files[block_id] >> min_premer_in_this_block;
		//second value in pair keeps track of the file from which the premer was drawn
		min_heap.push(pair<Premer<string>, uint64_t>(min_premer_in_this_block, block_id));
	}
    string output_path = (boost::format("%s/%s.out") % working_path % a_seq_local_path.filename().string()).str();
    BinaryDecoder bd;
	bd.set_k(k);
	ofstream out(output_path.c_str());
	while (min_heap.size() > 0) {
		pair<Premer<string>, uint64_t> min_premer_pair = min_heap.top();
		min_heap.pop();
		auto& a_premer = min_premer_pair.first;
		out << bd.decode(a_premer.m_key) << "\t" << a_premer.n_pos << "\t" << a_premer.m_val << "\n";
		Premer<string> cur_min_premer_in_this_block;
		if (input_files[min_premer_pair.second] >> cur_min_premer_in_this_block) {
			min_heap.push(pair<Premer<string>, uint64_t>(cur_min_premer_in_this_block, min_premer_pair.second));
		}
	}

	// close all input files
	for (uint64_t block_id = 0; block_id < max_block_index; ++block_id) {
		input_files[block_id].close();
	}
	//free input streams
	delete[] input_files;
	cout << checker;
	cout << "[MinimerSorter.sort_all_reads] wait until IO finishes\n";
}
bool MinimerSorter::process_a_fastq_line(ofstream& out, boost::mutex& a_mutex, const string& a_read, const string& a_qual) {
	out << a_read << "\n";
	return true;
}

bool MinimerSorter::process_a_fasta_line(ofstream& out, boost::mutex& a_mutex, const string& a_read) {
	out << a_read << "\n";
	return true;
}

void MinimerSorter::calculate_all_file_sizes() {
	TimeChecker checker;
	checker.setTarget("MinimerSorter.calculate_all_file_sizes");
	checker.start_without_output();
	file_sizes.resize(paths.size());
	file_types.resize(paths.size());

	vector<function<void()> > tasks;
	for (uint64_t a_path_id = 0; a_path_id < paths.size(); ++a_path_id) {
		tasks.push_back([&, a_path_id] {
			file_sizes[a_path_id] = IOUtils::get_file_size(paths[a_path_id]);
			file_types[a_path_id] = IOUtils::get_file_format(paths[a_path_id]);
		});
	}
	ParallelRunner::run_unbalanced_load(n_cores, tasks);
	for (uint64_t a_file_id = 0; a_file_id < file_sizes.size(); ++a_file_id) {
		total_file_size += file_sizes[a_file_id];
	}
	cout << (boost::format("[MinimerSorter.calculate_all_file_sizes] Total: %d bytes\n") % total_file_size).str();
	cout << checker;
}

void MinimerSorter::create_all_blocks() {
	TimeChecker checker;
	checker.setTarget("MinimerSorter.create_all_blocks");
	checker.start();
	n_blocks.resize(paths.size());
	skip_points.resize(paths.size());
	actual_chunk_sizes.resize(paths.size());

	int64_t max_using_memory_size = 4 * IOUtils::GIGA;
	double size_of_block = max_using_memory_size / n_cores;

	for (uint64_t a_path_id = 0; a_path_id < paths.size(); ++a_path_id) {
		n_blocks[a_path_id] = file_sizes[a_path_id] / (double) size_of_block;
		if(0 == n_blocks[a_path_id]) {
			n_blocks[a_path_id] = 1;
			vector<uint64_t>& a_set_of_skip_points = skip_points[a_path_id];
			a_set_of_skip_points.resize(2);
			a_set_of_skip_points[0] = 0;
			a_set_of_skip_points[1] = total_file_size;
			vector<uint64_t>& a_set_of_actual_chunk_sizes =
					actual_chunk_sizes[a_path_id];
			a_set_of_actual_chunk_sizes.resize(1);
			a_set_of_actual_chunk_sizes[0] = a_set_of_skip_points[1]
					- a_set_of_skip_points[0];
			continue;
		}
		if(n_blocks[a_path_id] < n_cores) {
			size_of_block = file_sizes[a_path_id] / (double) n_cores;
			n_blocks[a_path_id] = file_sizes[a_path_id] / (double) size_of_block;
		}

		vector<uint64_t>& a_set_of_skip_points = skip_points[a_path_id];
		a_set_of_skip_points.resize(n_blocks[a_path_id] + 1);
		if(file_types[a_path_id] == FASTQ) {
			IOUtils::find_fastq_skip_points(a_set_of_skip_points, paths[a_path_id],
				size_of_block, file_sizes[a_path_id], n_blocks[a_path_id],
				n_cores);
		} else if(file_types[a_path_id] == FASTA) {
			IOUtils::find_fasta_skip_points(a_set_of_skip_points, paths[a_path_id],
			size_of_block, file_sizes[a_path_id], n_blocks[a_path_id],
			n_cores);
		} else {
			cerr << "[MinimerSorter.create_all_blocks] Unsupported file format\n";
			exit(1);
		}
		vector<uint64_t>& a_set_of_actual_chunk_sizes =
				actual_chunk_sizes[a_path_id];
		a_set_of_actual_chunk_sizes.resize(n_blocks[a_path_id]);
		for (uint64_t skip_id = 0; skip_id < a_set_of_skip_points.size() - 1;
				++skip_id) {
			a_set_of_actual_chunk_sizes[skip_id] = a_set_of_skip_points[skip_id
					+ 1] - a_set_of_skip_points[skip_id];
		}
	}
	cout << checker;
}

void MinimerSorter::print_progress(boost::mutex& a_mutex, int64_t& n_processing) {
	boost::lock_guard<boost::mutex> lock(a_mutex);
	++n_processing;
	cout << "*";
	if (0 != n_processing && 0 == (n_processing % n_dividing_elements)) {
		cout << "\n";
	}
}

} /* namespace pillar */
