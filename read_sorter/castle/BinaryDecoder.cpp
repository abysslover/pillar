/*
 * BinaryDecoder.cpp
 *
 *  Created on: Sep 24, 2019
 *      Author: Eun-Cheon Lim @ Postech Plant Genomic Recombination Lab.
 *      Email : abysslover@gmail.com
 *    License : GPLv3
 */

#include "BinaryDecoder.hpp"

namespace castle {

char BinaryDecoder::decoding_map[] = { 'A', 'C', 'G', 'T' };

BinaryDecoder::BinaryDecoder() : k(10) {
}

BinaryDecoder::~BinaryDecoder() {
}
void BinaryDecoder::set_k(const uint64_t a_k) {
	k = a_k;
}
string BinaryDecoder::decode(const uint64_t a_val) {
	string result;
	uint32_t processed = 0;
	for (uint64_t shift_bits = TWO_BITS_INITIAL_SHIFT_BITS_LONG; processed < k && shift_bits >= 0; shift_bits -= TWO_BITS_UNIT_LONG) {
		 char c = decoding_map[(int) ((a_val & (TWO_BITS_MASK_LONG << shift_bits)) >> shift_bits)];
		 result += c;
		 ++processed;
	}
	return result;
}

} /* namespace castle */
