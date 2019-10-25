/*
 * SmallKmer.hpp
 *
 *  Created on: 2019. 10. 24.
 *      Author: pgr
 */

#ifndef PILLAR_SMALLKMER_HPP_
#define PILLAR_SMALLKMER_HPP_
#include <iostream>
#include <limits>
#include <bitset>

namespace pillar {
using namespace std;
class SmallKmer {
private:
	int32_t k;
	int32_t n_bits;
	int64_t bit_pos;
	int64_t term_bit_pos;
	bool is_max;
public:
	uint64_t m_bits;
	uint64_t m_reverse_bits;
	static uint8_t EncodingMap[256];
	static uint8_t ReverseEncodingMap[256];
	static unsigned char DecodingMap[6];
public:
	SmallKmer();
	~SmallKmer();
	void set_k(const int32_t a_k);
	void left_shift(const int32_t n_shift);
	void encode(const string& a_str, const uint64_t start_pos);
	void encode_a_letter(const char c, const char c_r);
	void encode_a_terminator();
	void encode_max();
	uint64_t get_minimizer();
	bool check_consecutive_letters(const uint64_t& bits);
	string decode() const;
	string decode_rev() const;
	void push_back(const char c);
	void shift_and_push_back(const char c);
	bool is_empty() const;
public:
	bool operator == (const SmallKmer &rhs) const {
		return this->m_bits == rhs.m_bits;
	}
	bool operator < (const SmallKmer &rhs) const {
		return this->m_bits < rhs.m_bits;
	}
	char operator [](int64_t i) const {
		int64_t local_bit_pos = i * 3;
		uint64_t cur_mask = 0b111;
		cur_mask <<= local_bit_pos;
		uint64_t b = m_bits & cur_mask;
		b >>= local_bit_pos;
		return DecodingMap[b];
	}
	char operator ()(int64_t i) const {
		int64_t local_bit_pos = i * 3;
		uint64_t cur_mask = 0b111;
		cur_mask <<= local_bit_pos;
		uint64_t b = m_reverse_bits & cur_mask;
		b >>= local_bit_pos;
		return DecodingMap[b];
	}
	size_t size() const {
		return k;
	}
};

} /* namespace pillar */

#endif /* PILLAR_SMALLKMER_HPP_ */
