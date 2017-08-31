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
#ifndef KMERCODE_HPP
#define KMERCODE_HPP 1

#include <cstdint>
#include <iostream>
#include "globals.hpp"
#include "murmur.hpp"
using namespace std;

typedef uint64_t block_t;
typedef uint32_t half_block_t;

static const int block_s = sizeof(uint64_t) * 8;

class KmerCode {

	private:
		block_t _code[Globals::NB_BLOCKS];
		

	public:
		KmerCode ();
		KmerCode (const unsigned int code);
		KmerCode (const KmerCode &code);

		KmerCode operator<<=(const unsigned int v);
		KmerCode operator>>=(const unsigned int v);
		KmerCode operator<<(const unsigned int v) const;
		KmerCode operator>>(const unsigned int v) const;

		KmerCode operator|=(const KmerCode &v);
		KmerCode operator&=(const KmerCode &v);
		KmerCode operator|(const KmerCode &v) const;
		KmerCode operator&(const KmerCode &v) const;

		unsigned int to_uint() const;

		uint32_t hash() const {
			//return murmurHash(_code, sizeof(block_t) * Globals::NB_BLOCKS, 0) + Globals::KEY_EMPTY + 1;
			return murmurHash(_code, sizeof(block_t) * Globals::NB_BLOCKS, 0);
		}

		friend bool operator== (const KmerCode &c1, const KmerCode &c2);
		friend bool operator!= (const KmerCode &c1, const KmerCode &c2);
		friend bool operator< (const KmerCode &c1, const KmerCode &c2);
		friend ostream& operator<< (ostream& output, const KmerCode &c);
};

namespace std {    
   template <>
   struct hash< KmerCode > {
       inline size_t operator()(const KmerCode v) const {
			return v.hash();
       }
   };
}

struct KmerCodeEqStr {
  bool operator()(const KmerCode &s1, const KmerCode &s2) const {
    return (s1 == s2);
  }
};


#endif
