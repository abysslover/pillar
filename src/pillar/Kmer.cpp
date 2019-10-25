/*
 * Kmer.cpp
 *
 *  Created on: 2019. 10. 22.
 *      Author: Eun-Cheon Lim @ Postech Plant Genomics Lab.
 */

#include "Kmer.hpp"

namespace pillar {
uint8_t Kmer::EncodingMap[] = {
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

uint8_t Kmer::ReverseEncodingMap[] = {
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
unsigned char Kmer::DecodingMap[] = { 'N', 'A', 'C', 'G', 'T', '$' };
Kmer::Kmer() : k(15), n_bits(45), bit_pos(0), term_bit_pos(0), is_max(false) {
}

Kmer::~Kmer() {
}

void Kmer::set_k(const int32_t a_k) {
	k = a_k;
	// for the '$' sign
	n_bits = 3 * (k + 1);
}
void Kmer::left_shift(const int32_t n_shift) {
	// remove terminate letter '$' from reverse complement
	m_reverse_bits[bit_pos] = 0;
	m_reverse_bits[bit_pos + 1] = 0;
	m_reverse_bits[bit_pos + 2] = 0;
	m_bits >>= (n_shift * 3);
	bit_pos -= 3;
	term_bit_pos = bit_pos;
	char c = '$';
	uint8_t b = EncodingMap[static_cast<int64_t>(c)];
	uint8_t first_b = b & 4;
	uint8_t center_b = b & 2;
	uint8_t last_b = b & 1;
	m_reverse_bits[term_bit_pos] = first_b;
	m_reverse_bits[term_bit_pos + 1] = center_b;
	m_reverse_bits[term_bit_pos + 2] = last_b;
	is_max = false;
}
void Kmer::encode(const string& a_str, const uint64_t start_pos) {
	m_bits.resize(n_bits);
	m_reverse_bits.resize(n_bits);
	uint64_t max_pos = min(a_str.size(), start_pos + k);
	bit_pos = 0;
	for(uint64_t b_id = start_pos; b_id < max_pos; ++b_id) {
		char c = a_str[b_id];
		uint8_t b = EncodingMap[static_cast<int64_t>(c)];
		uint8_t first_b = b & 4;
		uint8_t center_b = b & 2;
		uint8_t last_b = b & 1;
		m_bits[bit_pos] = first_b;
		m_bits[bit_pos + 1] = center_b;
		m_bits[bit_pos + 2] = last_b;

		char c_r = a_str[max_pos - b_id - 1];
		uint8_t r = ReverseEncodingMap[static_cast<int64_t>(c_r)];
		uint8_t first_r = r & 4;
		uint8_t center_r = r & 2;
		uint8_t last_r = r & 1;
		m_reverse_bits[bit_pos] = first_r;
		m_reverse_bits[bit_pos + 1] = center_r;
		m_reverse_bits[bit_pos + 2] = last_r;
		bit_pos += 3;
	}
	char c = '$';
	uint8_t b = EncodingMap[static_cast<int64_t>(c)];
	uint8_t first_b = b & 4;
	uint8_t center_b = b & 2;
	uint8_t last_b = b & 1;
	m_bits[bit_pos] = first_b;
	m_bits[bit_pos + 1] = center_b;
	m_bits[bit_pos + 2] = last_b;

	m_reverse_bits[bit_pos] = first_b;
	m_reverse_bits[bit_pos + 1] = center_b;
	m_reverse_bits[bit_pos + 2] = last_b;
	term_bit_pos = bit_pos;
	is_max = false;
}
void Kmer::encode_max() {
	m_bits.resize(n_bits);
	m_reverse_bits.resize(n_bits);
	m_bits.set();
	m_reverse_bits.set();
	is_max = true;
}
boost::dynamic_bitset<uint64_t> Kmer::get_minimizer() {
	bool has_consecutive_letters = check_consecutive_letters(m_bits);
	bool has_consecutive_reverse_letters = check_consecutive_letters(m_reverse_bits);
	if (has_consecutive_letters && has_consecutive_reverse_letters) {
		boost::dynamic_bitset<uint64_t> max_bits;
		max_bits.resize(n_bits);
		max_bits.set();
		return max_bits;
	}
	if(m_bits < m_reverse_bits) {
		return m_bits;
	} else if(m_bits > m_reverse_bits) {
		return m_reverse_bits;
	}
	return m_bits;
}

bool Kmer::check_consecutive_letters(const boost::dynamic_bitset<uint64_t>& bits) {
	int64_t b_id = 0;
	char prev_c = ' ';
	int64_t cur_len = 0;
	for(int64_t b_pos = 0; b_pos < n_bits; b_pos += 3, ++b_id) {
		uint8_t first_b = m_bits[b_pos];
		uint8_t center_b = m_bits[b_pos + 1];
		uint8_t last_b = m_bits[b_pos + 2];
		uint8_t temp_b = first_b << 2 | center_b << 1 | last_b;
		char c = DecodingMap[temp_b];
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
string Kmer::decode() const {
	string d;
	d.resize(k);
	int64_t b_id = 0;
	for(int64_t b_pos = 0; b_pos < n_bits; b_pos += 3, ++b_id) {
		uint8_t first_b = m_bits[b_pos];
		uint8_t center_b = m_bits[b_pos + 1];
		uint8_t last_b = m_bits[b_pos + 2];
		uint8_t temp_b = first_b << 2 | center_b << 1 | last_b;
		char c = DecodingMap[temp_b];
		if (c == '$') {
			break;
		}
		d[b_id] = c;
	}
	return d;
}

string Kmer::decode_rev() const {
	string d;
	d.resize(k);
	int64_t b_id = 0;
	for(int64_t b_pos = 0; b_pos < n_bits; b_pos += 3, ++b_id) {
		uint8_t first_b = m_reverse_bits[b_pos];
		uint8_t center_b = m_reverse_bits[b_pos + 1];
		uint8_t last_b = m_reverse_bits[b_pos + 2];
		uint8_t temp_b = first_b << 2 | center_b << 1 | last_b;
		char c = DecodingMap[temp_b];
		if (c == '$') {
			break;
		}
		d[b_id] = c;
	}
	return d;
}
void Kmer::push_back(const char c) {
	uint8_t b = EncodingMap[static_cast<int64_t>(c)];
	uint8_t first_b = b & 4;
	uint8_t center_b = b & 2;
	uint8_t last_b = b & 1;
//	cout << "Forward (Before): " << m_bits << "\n";
	m_bits[bit_pos] = first_b;
	m_bits[bit_pos + 1] = center_b;
	m_bits[bit_pos + 2] = last_b;

//	cout << "Reverse (Before): " << m_reverse_bits << "\n";

	uint8_t r = ReverseEncodingMap[static_cast<int64_t>(c)];

	uint8_t first_r = r & 4;
	uint8_t center_r = r & 2;
	uint8_t last_r = r & 1;

//	m_reverse_bits[bit_pos] = first_r;
//	m_reverse_bits[bit_pos + 1] = center_r;
//	m_reverse_bits[bit_pos + 2] = last_r;

	// replace the first character
	m_reverse_bits <<= 3;
	m_reverse_bits[0] = first_r;
	m_reverse_bits[1] = center_r;
	m_reverse_bits[2] = last_r;

	bit_pos += 3;
	char last_c = '$';
	b = EncodingMap[static_cast<int64_t>(last_c)];
	first_b = b & 4;
	center_b = b & 2;
	last_b = b & 1;
	m_bits[bit_pos] = first_b;
	m_bits[bit_pos + 1] = center_b;
	m_bits[bit_pos + 2] = last_b;
	m_reverse_bits[bit_pos] = first_b;
	m_reverse_bits[bit_pos + 1] = center_b;
	m_reverse_bits[bit_pos + 2] = last_b;
	term_bit_pos = bit_pos;
//	cout << "Forward (After) : " <<  m_bits << "\n";
//	cout << "Reverse (After) : " << m_reverse_bits << "\n";
	is_max = false;
}
void Kmer::shift_and_push_back(const char c) {
	left_shift(1);
	push_back(c);
}
bool Kmer::is_empty() const {
	return m_bits.none();
}
} /* namespace pillar */
