/*
 * Kmer.hpp
 *
 *  Created on: 2019. 10. 22.
 *      Author: Eun-Cheon Lim @ Postech Plant Genomics Lab.
 */

#ifndef PILLAR_KMER_HPP_
#define PILLAR_KMER_HPP_

#include <boost/dynamic_bitset.hpp>
#include <iostream>

namespace pillar {
using namespace std;
class Kmer {
private:
	int32_t k;
	int32_t n_bits;
	int64_t bit_pos;
	int64_t term_bit_pos;
	bool is_max;
public:
	boost::dynamic_bitset<uint64_t> m_bits;
	boost::dynamic_bitset<uint64_t> m_reverse_bits;
	static uint8_t EncodingMap[256];
	static uint8_t ReverseEncodingMap[256];
	static unsigned char DecodingMap[6];
public:
	Kmer();
	~Kmer();
	void set_k(const int32_t a_k);
	void left_shift(const int32_t n_shift);
	void encode(const string& a_str, const uint64_t start_pos);
	void encode_max();
	boost::dynamic_bitset<uint64_t> get_minimizer();
	bool check_consecutive_letters(const boost::dynamic_bitset<uint64_t>& bits);
	string decode() const;
	string decode_rev() const;
	void push_back(const char c);
	void shift_and_push_back(const char c);
	bool is_empty() const;
public:
	bool operator == (const Kmer &rhs) const {
		return this->m_bits == rhs.m_bits;
	}
	bool operator < (const Kmer &rhs) const {
		return this->m_bits < rhs.m_bits;
	}
	char operator [](int64_t i) const {
		int64_t local_bit_pos = i * 3;
		uint8_t first_b = m_bits[local_bit_pos];
		uint8_t center_b = m_bits[local_bit_pos + 1];
		uint8_t last_b = m_bits[local_bit_pos + 2];
		uint8_t temp_b = first_b << 2 | center_b << 1 | last_b;
		return DecodingMap[temp_b];
	}
	char operator ()(int64_t i) const {
		int64_t local_bit_pos = i * 3;
		uint8_t first_b = m_reverse_bits[local_bit_pos];
		uint8_t center_b = m_reverse_bits[local_bit_pos + 1];
		uint8_t last_b = m_reverse_bits[local_bit_pos + 2];
		uint8_t temp_b = first_b << 2 | center_b << 1 | last_b;
		return DecodingMap[temp_b];
	}
	size_t size() const {
		return k;
	}
};

} /* namespace pillar */

#endif /* PILLAR_KMER_HPP_ */
