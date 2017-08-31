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
#ifndef HASHES_HPP
#define HASHES_HPP 1

#include "globals.hpp"
#include "kmerCode.hpp"

#ifdef HASH_SLOW
#include "sparsehash/sparse_hash_map"
using google::sparse_hash_map;
typedef sparse_hash_map<KmerCode, KmerNb, hash<KmerCode>, KmerCodeEqStr> _SparseHash;
class SparseHash: public _SparseHash {
	public:
		SparseHash() {
			set_deleted_key(Globals::KEY_DELETED);
		}
};
#endif
#ifdef HASH_MID
#include <unordered_map>
using namespace std;
typedef unordered_map<KmerCode, KmerNb, hash<KmerCode>, KmerCodeEqStr> SimpleHash;
#endif
#ifdef HASH_FAST
#include "sparsehash/dense_hash_map"
using google::dense_hash_map;
typedef dense_hash_map<KmerCode, KmerNb, hash<KmerCode>, KmerCodeEqStr> _DenseHash;
class DenseHash: public _DenseHash {
	public:
		DenseHash() {
			set_empty_key(Globals::KEY_EMPTY);
			set_deleted_key(Globals::KEY_DELETED);
		}
};
#endif

#endif
