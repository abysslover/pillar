/*
 * SplitReader.cpp
 *
 *  Created on: 2019. 10. 22.
 *      Author: Eun-Cheon Lim @ Postech Plant Genomics Lab.
 */

#include "SplitReader.hpp"

namespace pillar {

SplitReader::SplitReader() : n_dividing_elements(100), total_file_size(0), window_size(101), minimizer_size(7){
	TimeChecker checker;
	n_cores = checker.get_number_of_cores();
}

SplitReader::~SplitReader() {
}
void SplitReader::split_a_read(const string& a_read, const string& a_qual) {
	cout << a_read << "\n";
	int64_t max_id = a_read.size();
	max_id -= minimizer_size;
	SmallKmer kmer;
	kmer.set_k(minimizer_size);
	vector<SmallKmer> kmers;
	for(int64_t b_id = 0; b_id <= max_id; ++b_id) {
		if(kmer.is_empty()) {
			kmer.encode(a_read, b_id);
		} else {
			int64_t cur_id = b_id + minimizer_size - 1;
			char c = a_read[cur_id];
			kmer.shift_and_push_back(c);
		}
		kmers.push_back(kmer);
	}
	max_id = a_read.size();
	max_id -= window_size;
	SmallKmer prev_minimizer;
	prev_minimizer.set_k(minimizer_size);
	prev_minimizer.encode_max();
	int64_t prev_id = 0;
	int64_t state = 0;
	for(int64_t b_id = 0; b_id <= max_id; ++b_id) {
		auto a_minimizer = get_minimizer(kmers, b_id);
		switch(state) {
		case 0:
			prev_minimizer = a_minimizer;
			state = 1;
			continue;
		case 1:
			if(prev_minimizer.get_minimizer() != a_minimizer.get_minimizer()) {
//				if(0 == b_id) {
//	//				auto cur_str = a_read.substr(b_id, window_size);
//	//				cout << "NORMAL: " << cur_str << "\t" << cur_str.size() << "\t" << a_minimizer.decode() << "\t" << a_minimizer.decode_rev() << "\n";
//				} else {
//					cout << "BREAK : " << b_id << "\t" << a_minimizer.decode() << "\t" << a_minimizer.decode_rev() << "\n";
					int64_t cur_len = b_id - prev_id + window_size - 1;
					cout << "TRANSI: " << b_id << "\t" << a_read.substr(prev_id, cur_len) << "\t" << cur_len << "\t" << prev_minimizer.decode() << "\t" << a_minimizer.decode_rev() << "\n";
					b_id = prev_id + cur_len;
					state = 0;
//				}
				prev_id = b_id;
				prev_minimizer = a_minimizer;
			}
//			else {
//				auto cur_str = a_read.substr(b_id, window_size);
//				cout << "NORMAL: " << cur_str << "\t" << cur_str.size() << "\t" << a_minimizer.decode() << "\t" << a_minimizer.decode_rev() << "\n";
//			}
		}
	}
}

void SplitReader::split_a_read_alt(const string& a_read, const string& a_qual) {
	cout << a_read << "\n";
	int64_t max_id = a_read.size();
	max_id -= minimizer_size;
	Kmer kmer;
	kmer.set_k(minimizer_size);
	vector<Kmer> kmers;
	for(int64_t b_id = 0; b_id <= max_id; ++b_id) {
		if(kmer.is_empty()) {
			kmer.encode(a_read, b_id);
		} else {
			int64_t cur_id = b_id + minimizer_size - 1;
			char c = a_read[cur_id];
			kmer.shift_and_push_back(c);
		}
		kmers.push_back(kmer);
	}
	max_id = a_read.size();
	max_id -= window_size;
	Kmer prev_minimizer;
	prev_minimizer.set_k(minimizer_size);
	prev_minimizer.encode_max();
	int64_t prev_id = 0;
	int64_t state = 0;
	for(int64_t b_id = 0; b_id <= max_id; ++b_id) {
		auto a_minimizer = get_minimizer_alt(kmers, b_id);
		switch(state) {
		case 0:
			prev_minimizer = a_minimizer;
			state = 1;
			continue;
		case 1:
			if(prev_minimizer.get_minimizer() != a_minimizer.get_minimizer()) {
				if(0 == b_id) {
	//				auto cur_str = a_read.substr(b_id, window_size);
	//				cout << "NORMAL: " << cur_str << "\t" << cur_str.size() << "\t" << a_minimizer.decode() << "\t" << a_minimizer.decode_rev() << "\n";
				} else {
//					cout << "BREAK : " << b_id << "\t" << a_minimizer.decode() << "\t" << a_minimizer.decode_rev() << "\n";
					int64_t cur_len = b_id - prev_id + window_size - 1;
					cout << "TRANSI: " << b_id << "\t" << a_read.substr(prev_id, cur_len) << "\t" << cur_len << "\t" << prev_minimizer.decode() << "\t" << a_minimizer.decode_rev() << "\n";
					b_id = prev_id + cur_len;
					state = 0;
				}
				prev_id = b_id;
				prev_minimizer = a_minimizer;
			}
//			else {
//				auto cur_str = a_read.substr(b_id, window_size);
//				cout << "NORMAL: " << cur_str << "\t" << cur_str.size() << "\t" << a_minimizer.decode() << "\t" << a_minimizer.decode_rev() << "\n";
//			}
		}
	}
}
SmallKmer SplitReader::get_single_minimizer(vector<SmallKmer>& kmers, const uint64_t start_id) {
	uint64_t max_id = min(kmers.size(), start_id + window_size);
	SmallKmer a_minimizer;
	a_minimizer.encode_max();

	for(uint64_t k_id = start_id; k_id < max_id; ++k_id) {
		auto& a_kmer = kmers[k_id];
		auto cur_minimizer = a_kmer.get_minimizer();
		auto min_minimizer = a_minimizer.get_minimizer();
		if (cur_minimizer < min_minimizer) {
			a_minimizer.m_bits = cur_minimizer;
		}
	}
	return a_minimizer;
}
SmallKmer SplitReader::get_minimizer(vector<SmallKmer>& kmers, const uint64_t start_id) {
	int64_t n_kmer_based_max_id = kmers.size();
	n_kmer_based_max_id -= window_size;
	int64_t sliding_window_based_max_id = start_id + window_size;
	uint64_t max_id = min(n_kmer_based_max_id, sliding_window_based_max_id);
	SmallKmer a_minimizer;
	a_minimizer.encode_max();

	for(uint64_t k_id = start_id; k_id < max_id; ++k_id) {
		auto& a_kmer = kmers[k_id];
		auto cur_minimizer = a_kmer.get_minimizer();
		auto min_minimizer = a_minimizer.get_minimizer();
		if (cur_minimizer < min_minimizer) {
			a_minimizer = a_kmer;
		}
	}
	return a_minimizer;
}
Kmer SplitReader::get_minimizer_alt(vector<Kmer>& kmers, const uint64_t start_id) {
	int64_t n_kmer_based_max_id = kmers.size();
	n_kmer_based_max_id -= window_size;
	int64_t sliding_window_based_max_id = start_id + window_size;
	uint64_t max_id = min(n_kmer_based_max_id, sliding_window_based_max_id);
	Kmer a_minimizer;
	a_minimizer.encode_max();

	for(uint64_t k_id = start_id; k_id < max_id; ++k_id) {
		auto& a_kmer = kmers[k_id];
		auto cur_minimizer = a_kmer.get_minimizer();
		auto min_minimizer = a_minimizer.get_minimizer();
		if (cur_minimizer < min_minimizer) {
			a_minimizer = a_kmer;
		}
	}
	return a_minimizer;
}
void SplitReader::add_read_path(const string& a_path) {
	if (!bfs::exists(a_path)) {
		return;
	}
	paths.push_back(a_path);
}
void SplitReader::set_working_path(const string& a_path) {
	working_path = a_path;
	if (!bfs::exists(a_path)) {
		bfs::create_directories(a_path);
	}
}
void SplitReader::prepare_parallel_processing() {
	calculate_all_file_sizes();
	create_all_blocks();
}
void SplitReader::calculate_all_file_sizes() {
	TimeChecker checker;
	checker.setTarget("SplitReader.calculate_all_file_sizes");
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
	cout << (boost::format("[SplitReader.calculate_all_file_sizes] Total: %d bytes\n") % total_file_size).str();
	cout << checker;
}

void SplitReader::create_all_blocks() {
	TimeChecker checker;
	checker.setTarget("SplitReader.create_all_blocks");
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
			cerr << "[SplitReader.create_all_blocks] Unsupported file format\n";
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
void SplitReader::create_splits() {
	for(uint64_t a_seq_id = 0; a_seq_id < paths.size(); ++ a_seq_id) {
		find_a_split(a_seq_id);
	}
}

void SplitReader::find_a_split(const int64_t the_seq_id) {
	TimeChecker checker;
	checker.setTarget("SplitReader.find_a_split");
	checker.start_without_output();
	vector<function<void()> > tasks;
	string a_seq_path = paths[the_seq_id];
	const uint64_t max_block_index = actual_chunk_sizes[the_seq_id].size();
	cout << (boost::format("[SplitReader.find_a_split] # blocks: %d\n") % max_block_index).str();
	boost::mutex print_mutex;
	boost::filesystem::path a_seq_local_path(a_seq_path);
	vector<uint64_t> out_n_succeeded_count(max_block_index + 1);
	vector<uint64_t> out_n_processed_count(max_block_index + 1);

	set<uint64_t> debug_block;
	debug_block.insert(0);
	int64_t n_processing = 0;
	for (uint64_t block_id = 0; block_id < max_block_index; ++block_id) {
		tasks.push_back([&, &n_processing, the_seq_id, block_id] {
			if(debug_block.end() == debug_block.find(block_id)) {
//				print_progress(print_mutex, n_processing);
				return;
			}
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
						split_a_read(a_read, a_qual);

						if (a_local_processed_bytes > IOUtils::PARTIAL_READ_BUFFER_SIZE) {
							s_read.str(string());
							s_qual.str(string());
							a_local_processed_bytes = 0;
						}

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
					split_a_read(a_read, a_qual);
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
						string a_qual;
						split_a_read(a_read, a_qual);

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
					split_a_read(a_read, a_qual);
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
			print_progress(print_mutex, n_processing);
		});
	}
	ParallelRunner::run_unbalanced_load(n_cores, tasks);
}

void SplitReader::find_all_minimizers() {
	const function<void(set<SmallKmer>&, boost::mutex&, string&, string&)> a_callback = [&](set<SmallKmer>& out, boost::mutex& print_mutex, string& a_seq, string& a_qual){
		find_binning_minimizers(out, a_seq);
	};
	for(uint64_t a_seq_id = 0; a_seq_id < paths.size(); ++ a_seq_id) {
		set<SmallKmer> final_results;
		auto results = process_a_block_task(a_seq_id, "find_minimizers_in_block", a_callback);
		for(auto& result_set: results) {
			final_results.insert(result_set.begin(), result_set.end());
		}
		cout << (boost::format("[SplitReader::find_all_minimizers] %s: %d\n") % paths[a_seq_id] % final_results.size()).str();
	}
}

void SplitReader::find_binning_minimizers(set<SmallKmer>& result_minimizers, const string& a_read) {
//	cout << a_read << "\n";
	int64_t max_id = a_read.size();
	max_id -= minimizer_size;
	SmallKmer kmer;
	kmer.set_k(minimizer_size);
//	set<Kmer> result_minimizers;
	vector<SmallKmer> kmers;
	for(int64_t b_id = 0; b_id <= max_id; ++b_id) {
		if(kmer.is_empty()) {
			kmer.encode(a_read, b_id);
		} else {
			int64_t cur_id = b_id + minimizer_size - 1;
			char c = a_read[cur_id];
			kmer.shift_and_push_back(c);
		}
		kmers.push_back(kmer);
	}
	max_id = a_read.size();
	max_id -= window_size;
	SmallKmer prev_minimizer;
	prev_minimizer.set_k(minimizer_size);
	prev_minimizer.encode_max();
	int64_t prev_id = 0;
	int64_t state = 0;
	for(int64_t b_id = 0; b_id <= max_id; ++b_id) {
		auto a_minimizer = get_single_minimizer(kmers, b_id);
		switch(state) {
		case 0:
			prev_minimizer = a_minimizer;
			state = 1;
			continue;
		case 1:
			if(prev_minimizer.m_bits != a_minimizer.m_bits) {
				if(0 == b_id) {
	//				auto cur_str = a_read.substr(b_id, window_size);
	//				cout << "NORMAL: " << cur_str << "\t" << cur_str.size() << "\t" << a_minimizer.decode() << "\t" << a_minimizer.decode_rev() << "\n";
				} else {
//					cout << "BREAK : " << b_id << "\t" << a_minimizer.decode() << "\t" << a_minimizer.decode_rev() << "\n";
					int64_t cur_len = b_id - prev_id + window_size - 1;
//					cout << "TRANSI: " << b_id << "\t" << a_read.substr(prev_id, cur_len) << "\t" << cur_len << "\t" << a_minimizer.decode() << "\t" << a_minimizer.decode_rev() << "\n";
					b_id = prev_id + cur_len;
					result_minimizers.insert(a_minimizer);
					state = 0;
				}
				prev_id = b_id;
				prev_minimizer = a_minimizer;
			}
//			else {
//				auto cur_str = a_read.substr(b_id, window_size);
//				cout << "NORMAL: " << cur_str << "\t" << cur_str.size() << "\t" << a_minimizer.decode() << "\t" << a_minimizer.decode_rev() << "\n";
//			}
		}
	}
//	return result_minimizers;
}

void SplitReader::print_progress(boost::mutex& a_mutex, int64_t& n_processing) {
	boost::lock_guard<boost::mutex> lock(a_mutex);
	++n_processing;
	cout << "*";
	if (0 != n_processing && 0 == (n_processing % n_dividing_elements)) {
		cout << "\n";
	}
}
} /* namespace pillar */
