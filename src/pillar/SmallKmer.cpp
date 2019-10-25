/*
 * SmallKmer.cpp
 *
 *  Created on: 2019. 10. 24.
 *      Author: pgr
 */

#include "SmallKmer.hpp"

namespace pillar {

uint8_t SmallKmer::EncodingMap[] = {
	// 0 - 15
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	// 16 - 31
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	// 32 - 47
	0, 0, 0, 0, 5 /* '$' */, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	// 48 - 63
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	// 64 - 79
	0, 1 /* 'A' */, 0, 2 /* 'C' */, 0, 0, 0, 3 /* 'G' */, 0, 0, 0, 0, 0, 0, 0 /* 'N' */, 0,
	// 80 - 95
	0, 0, 0, 0, 4 /* 'T' */, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	// 96 - 111
	0, 1 /* 'a' */, 0, 2 /* 'c' */, 0, 0, 0, 3 /* 'g' */, 0, 0, 0, 0, 0, 0, 0 /* 'n' */, 0,
	// 112 - 127
	0, 0, 0, 0, 4 /* 't' */, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	// 128 - 143
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	// 144 - 159
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	// 160 - 175
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	// 176 - 191
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	// 192 - 207
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	// 208 - 223
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	// 224 - 239
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	// 240 - 255
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
};

uint8_t SmallKmer::ReverseEncodingMap[] = {
	// 0 - 15
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	// 16 - 31
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	// 32 - 47
	0, 0, 0, 0, 5 /* '$' */, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	// 48 - 63
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	// 64 - 79
	0, 4 /* 'A' */, 0, 3 /* 'C' */, 0, 0, 0, 2 /* 'G' */, 0, 0, 0, 0, 0, 0, 0 /* 'N' */, 0,
	// 80 - 95
	0, 0, 0, 0, 1 /* 'T' */, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	// 96 - 111
	0, 4 /* 'a' */, 0, 3 /* 'c' */, 0, 0, 0, 2 /* 'g' */, 0, 0, 0, 0, 0, 0, 0 /* 'n' */, 0,
	// 112 - 127
	0, 0, 0, 0, 1 /* 't' */, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	// 128 - 143
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	// 144 - 159
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	// 160 - 175
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	// 176 - 191
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	// 192 - 207
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	// 208 - 223
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	// 224 - 239
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	// 240 - 255
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
};
unsigned char SmallKmer::DecodingMap[] = { 'N', 'A', 'C', 'G', 'T', '$' };
SmallKmer::SmallKmer() : k(15), n_bits(45), bit_pos(0), term_bit_pos(0), is_max(false), m_bits(0), m_reverse_bits(0) {
}

SmallKmer::~SmallKmer() {
}
void SmallKmer::set_k(const int32_t a_k) {
	k = a_k;
	if (k > 21) {
		k = 21;
	}
	// for the '$' sign
	n_bits = 3 * (k + 1);
}
void SmallKmer::left_shift(const int32_t n_shift) {
//	cout << "Forward(Before): " << bitset<64>(m_bits) << "\n";
//	cout << "Reverse(Before): " << bitset<64>(m_reverse_bits) << "\n";
	// remove terminate letter '$' from reverse complement
	uint64_t cur_mask = 0b111;
	cur_mask <<= bit_pos;
	//e.g. inv_mask: 1111000111
	uint64_t inv_mask = ~cur_mask;
	m_reverse_bits &= inv_mask;

	bit_pos -= 3;
	term_bit_pos = bit_pos;
	m_bits >>= (n_shift * 3);
	char c = '$';
	uint64_t b = EncodingMap[static_cast<int64_t>(c)];
	//e.g. cur_mask: 0000111000
	cur_mask = 0b111;
	cur_mask <<= bit_pos;
	//e.g. inv_mask: 1111000111
	inv_mask = ~cur_mask;
	m_reverse_bits &= inv_mask;
	b <<= bit_pos;
	m_reverse_bits |= b;
//	cout << "Forward(After) : " << bitset<64>(m_bits) << "\n";
//	cout << "Reverse(After) : " << bitset<64>(m_reverse_bits) << "\n";

	is_max = false;
}
void SmallKmer::encode(const string& a_str, const uint64_t start_pos) {
	uint64_t max_pos = min(a_str.size(), start_pos + k);
	bit_pos = 0;
	for(uint64_t b_id = start_pos; b_id < max_pos; ++b_id) {
		char c = a_str[b_id];
		char c_r = a_str[max_pos - b_id - 1];
		encode_a_letter(c, c_r);
	}
	encode_a_terminator();
	is_max = false;
}
void SmallKmer::encode_a_letter(const char c, const char c_r) {
	uint64_t b = EncodingMap[static_cast<int64_t>(c)];
//	cout << "Forward(Before): " << bitset<64>(m_bits) << "\n";
	//e.g. cur_mask: 0000111000
	uint64_t cur_mask = 0b111;
	cur_mask <<= bit_pos;
	//e.g. inv_mask: 1111000111
	uint64_t inv_mask = ~cur_mask;

	m_bits &= inv_mask;
	b <<= bit_pos;
	m_bits |= b;
//	cout << "Forward(After) : " << bitset<64>(m_bits) << "\n";
//	cout << "Reverse(Before): " << bitset<64>(m_reverse_bits) << "\n";

	uint64_t r = ReverseEncodingMap[static_cast<int64_t>(c_r)];
	m_reverse_bits &= inv_mask;
	r <<= bit_pos;
	m_reverse_bits |= r;
//	cout << "Reverse(After) : " << bitset<64>(m_reverse_bits) << "\n";
	bit_pos += 3;
}
void SmallKmer::encode_a_terminator() {
	const char c = '$';
	uint64_t b = EncodingMap[static_cast<int64_t>(c)];
	//e.g. cur_mask: 0000111000
	uint64_t cur_mask = 0b111;
	cur_mask <<= bit_pos;
	//e.g. inv_mask: 1111000111
	uint64_t inv_mask = ~cur_mask;

	m_bits &= inv_mask;
	b <<= bit_pos;
	m_bits |= b;

	m_reverse_bits &= inv_mask;
	m_reverse_bits |= b;
	term_bit_pos = bit_pos;
}
void SmallKmer::encode_max() {
	m_reverse_bits = m_bits = UINT64_MAX;
	is_max = true;
}
uint64_t SmallKmer::get_minimizer() {
	bool has_consecutive_letters = check_consecutive_letters(m_bits);
	bool has_consecutive_reverse_letters = check_consecutive_letters(m_reverse_bits);
	if (has_consecutive_letters && has_consecutive_reverse_letters) {
		uint64_t max_bits = UINT64_MAX;
		return max_bits;
	}
	if(m_bits < m_reverse_bits) {
		return m_bits;
	} else if(m_bits > m_reverse_bits) {
		return m_reverse_bits;
	}
	return m_bits;
}

bool SmallKmer::check_consecutive_letters(const uint64_t& bits) {
	int64_t b_id = 0;
	char prev_c = ' ';
	int64_t cur_len = 0;
	for(int64_t b_pos = 0; b_pos < n_bits; b_pos += 3, ++b_id) {
		uint64_t cur_mask = 0b111;
		cur_mask <<= b_pos;
		uint64_t b = bits & cur_mask;
		b >>= b_pos;
		char c = DecodingMap[b];
		if (c == prev_c) {
			++cur_len;
			if (cur_len > 1) {
				break;
			}
		}
		prev_c = c;
		if (c == '$') {
			break;
		}
	}
	return cur_len > 1;
}
string SmallKmer::decode() const {
	string d;
	d.resize(k);
	int64_t b_id = 0;
	for(int64_t b_pos = 0; b_pos < n_bits; b_pos += 3, ++b_id) {
		uint64_t cur_mask = 0b111;
		cur_mask <<= b_pos;
		uint64_t b = m_bits & cur_mask;
		b >>= b_pos;
		char c = DecodingMap[b];
		if (c == '$') {
			break;
		}
		d[b_id] = c;
	}
	return d;
}

string SmallKmer::decode_rev() const {
	string d;
	d.resize(k);
	int64_t b_id = 0;
	for(int64_t b_pos = 0; b_pos < n_bits; b_pos += 3, ++b_id) {
		uint64_t cur_mask = 0b111;
		cur_mask <<= b_pos;
		uint64_t b = m_reverse_bits & cur_mask;
		b >>= b_pos;
		char c = DecodingMap[b];
		if (c == '$') {
			break;
		}
		d[b_id] = c;
	}
	return d;
}
void SmallKmer::push_back(const char c) {
//	cout << "Forward(Before): " << bitset<64>(m_bits) << "\n";
//	cout << "Reverse(Before): " << bitset<64>(m_reverse_bits) << "\n";
	uint64_t b = EncodingMap[static_cast<int64_t>(c)];
	uint64_t cur_mask = 0b111;
	cur_mask <<= bit_pos;
	//e.g. inv_mask: 1111000111
	uint64_t inv_mask = ~cur_mask;

	m_bits &= inv_mask;
	b <<= bit_pos;
	m_bits |= b;

	uint64_t r = ReverseEncodingMap[static_cast<int64_t>(c)];
	m_reverse_bits <<= 3;
	m_reverse_bits |= r;
	bit_pos += 3;
	encode_a_terminator();
//	cout << "Forward(After) : " << bitset<64>(m_bits) << "\n";
//	cout << "Reverse(After) : " << bitset<64>(m_reverse_bits) << "\n";
	is_max = false;
}
void SmallKmer::shift_and_push_back(const char c) {
	left_shift(1);
	push_back(c);
}
bool SmallKmer::is_empty() const {
	return 0 == m_bits;
}
} /* namespace pillar */
