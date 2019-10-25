/*
 * BinaryDecoder.hpp
 *
 *  Created on: Sep 24, 2019
 *      Author: Eun-Cheon Lim @ Postech Plant Genomic Recombination Lab.
 *      Email : abysslover@gmail.com
 *    License : GPLv3
 */

#ifndef CASTLE_BINARYDECODER_HPP_
#define CASTLE_BINARYDECODER_HPP_
#include <sstream>
#include <string>
#include <tuple>
#include <boost/format.hpp>
#include <boost/cstdint.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/multiprecision/cpp_int.hpp>

namespace castle {
using namespace std;
class BinaryDecoder {
private:
	uint64_t k;
public:
	static const uint64_t TWO_BITS_UNIT_LONG = 2;
	static const uint64_t TWO_BITS_INITIAL_SHIFT_BITS_LONG = 62;
	static const uint64_t TWO_BITS_MASK_LONG = 0b11;
	static char decoding_map[4];
public:

	BinaryDecoder();
	~BinaryDecoder();
	void set_k(const uint64_t a_k);
	string decode(const uint64_t a_val);
};

} /* namespace castle */

#endif /* CASTLE_BINARYDECODER_HPP_ */
