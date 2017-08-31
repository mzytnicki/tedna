/**
Copyright (C) 2013 INRA-URGI
This file is part of TEDNA, a short reads transposable elements assembler
TEDNA is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
GNU Affero General Public License for more details.
See the GNU Affero General Public License for more details.
You should have received a copy of the GNU Affero General Public License
along with this program.
**/
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "kmerCode.hpp"

KmerCode::KmerCode() {
	for (int i = 0; i < Globals::NB_BLOCKS; i++) {
		_code[i] = -1;
	}
}

KmerCode::KmerCode(const unsigned int code) {
	_code[0] = code;
	for (int i = 1; i < Globals::NB_BLOCKS; i++) {
		_code[i] = 0;
	}
}

KmerCode::KmerCode(const KmerCode &code) {
	for (int i = 0; i < Globals::NB_BLOCKS; i++) {
		_code[i] = code._code[i];
	}
}

KmerCode KmerCode::operator|=(const KmerCode &v) {
	for (int i = 0; i < Globals::NB_BLOCKS; i++) {
		_code[i] |= v._code[i];
	}
	return *this;
}

KmerCode KmerCode::operator&=(const KmerCode &v) {
	for (int i = 0; i < Globals::NB_BLOCKS; i++) {
		_code[i] &= v._code[i];
	}
	return *this;
}

KmerCode KmerCode::operator|(const KmerCode &v) const {
	return KmerCode(*this) |= v;
}

KmerCode KmerCode::operator&(const KmerCode &v) const {
	return KmerCode(*this) &= v;
}

KmerCode KmerCode::operator<<=(const unsigned int v) {
	block_t newBlock[Globals::NB_BLOCKS];
	for (int i = 0; i < Globals::NB_BLOCKS; i++) {
		newBlock[i] = 0;
	}
	for (long i = Globals::NB_BLOCKS-1; i >= v/block_s; i--) {
		newBlock[i] = _code[i - v/block_s] << (v % block_s);
		if ((i > v/block_s) && (v % block_s != 0)) {
			newBlock[i] |= _code[i - v/block_s - 1] >> (block_s - v % block_s);
		}
	}
	for (int i = 0; i < Globals::NB_BLOCKS; i++) {
		_code[i] = newBlock[i];
	}
	return *this;
}

KmerCode KmerCode::operator>>=(const unsigned int v) {
	block_t newBlock[Globals::NB_BLOCKS];
	for (int i = 0; i < Globals::NB_BLOCKS; i++) {
		newBlock[i] = 0;
	}
	for (unsigned int i = 0; i < Globals::NB_BLOCKS - v/block_s; i++) {
		newBlock[i] = _code[i + v/block_s] >> (v % block_s);
		if ((i + v/block_s < Globals::NB_BLOCKS-1) && (v % block_s != 0)) {
			newBlock[i] |= _code[i + v/block_s + 1] << (block_s - v % block_s);
		}
	}
	for (int i = 0; i < Globals::NB_BLOCKS; i++) {
		_code[i] = newBlock[i];
	}
	return *this;
}

KmerCode KmerCode::operator<<(const unsigned int v) const {
	return KmerCode(*this) <<= v;
}

KmerCode KmerCode::operator>>(const unsigned int v) const {
	return KmerCode(*this) >>= v;
}

/*
template <int _kmerCodeSize>
KmerCode KmerCode::operator+=(const unsigned int v) {
	unsigned int carry = v;
	block_t     top   = -1;
	for (int i = 0; i < _kmerCodeSize; i++) {
		_code[i] += carry;
		carry     = (top - _code[i] > carry)? 1: 0;
	}
	return *this;
}

template <int _kmerCodeSize>
KmerCode KmerCode::operator*=(const unsigned int v) {
	half_block_t carry = 0;
	half_block_t vl    = 0, vr = v;
	for (int i = 0; i < _kmerCodeSize; i++) {
		half_block_t cl = _code[i] >> 32, cr = _code[i];
		block_t rl = cl * vl, rm = cr * vl + cl * vr, rr = cr * vr + carry;
		_code[i] = rm << 32 + rr;
		carry    = rl + rr >> 32;
	}
	return *this;
}

template <int _kmerCodeSize>
KmerCode KmerCode::operator/=(const unsigned int v) {
	block_t rest = 0;
	for (int i = _kmerCodeSize-1; i >= 0; i--) {
		block_t c = _code[i] + (((rest << 32) / v) << 32);
		_code[i]   = c / v;
		rest       = c % v;
	}
	return *this;
}

template <int _kmerCodeSize>
KmerCode KmerCode::operator+(const unsigned int v) {
	return KmerCode(*this) += v;
}

template <int _kmerCodeSize>
KmerCode KmerCode::operator*(const unsigned int v) {
	return KmerCode(*this) *= v;
}

template <int _kmerCodeSize>
KmerCode KmerCode::operator/(const unsigned int v) {
	return KmerCode(*this) /= v;
}

template <int _kmerCodeSize>
unsigned int KmerCode::operator%(const unsigned int v) {
	block_t top        = -1;
	unsigned int base   = ((top % v) + 1) % v;
	unsigned int modulo = 0;
	unsigned int c      = 1;
	for (int i = 1; i < _kmerCodeSize; i++) {
		modulo += ((_code[i] % v)  * c) % v;
		c      *= (c * base) % v;
		if (c == 0) {
			return modulo;
		}
	}
	return modulo;
}

template <int _kmerCodeSize>
KmerCode pow(const KmerCode &c, const int i) {
	int s      = 1;
	int d      = 0;
	KmerCode r = 1;
	KmerCode p = 1;
	while (d < i) {
		if (i & s == 1) {
			r *= p;
			d += s;
		}
		p *= c;
		s <<= 1;
	}
	return (r);
}

*/

unsigned int KmerCode::to_uint() const {
	return static_cast<unsigned int>(_code[0]);
}

/*
inline size_t KmerCode::hash() const {
	size_t sum = 0;
	for (int i = 0; i < Globals::NB_BLOCKS; i++) {
		sum += _code[i];
	}
	return sum;
}
*/

bool operator==(const KmerCode &c1, const KmerCode &c2) {
	for (int i = 0; i < Globals::NB_BLOCKS; i++) {
		if (c1._code[i] != c2._code[i]) {
			return false;
		}
	}
	return true;
}

bool operator!=(const KmerCode &c1, const KmerCode &c2) {
	return ! (c1 == c2);
}

bool operator<(const KmerCode &c1, const KmerCode &c2) {
	for (int i = Globals::NB_BLOCKS-1; i >= 0; i--) {
		if (c1._code[i] < c2._code[i]) {
			return true;
		}
		if (c1._code[i] > c2._code[i]) {
			return false;
		}
	}
	return false;
}

ostream& operator<<(ostream& output, const KmerCode &c) {
	output << c._code[Globals::NB_BLOCKS-1];
	for (int i = Globals::NB_BLOCKS-2; i >= 0; i--) {
		output << "," << c._code[i];
	}
	return output;
}

/*
namespace std {    
   template <>
   struct hash< KmerCode > {
       inline size_t operator()(const KmerCode v) const {
           return v.hash();
       }
   };
}
*/
