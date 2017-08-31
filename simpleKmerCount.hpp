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
#ifndef SIMPLE_KMER_COUNT_HPP
#define SIMPLE_KMER_COUNT_HPP 1

#include <unordered_map>
#include <vector>
#include <mutex>
#include "globals.hpp"
#include "hashes.hpp"
#include "kmer.hpp"
using namespace std;

class SimpleKmerCount {

    protected:
#ifdef HASH_MID
		SimpleHash        _map;
#endif
#ifdef HASH_SLOW
		SparseHash        _map;
#endif
#ifdef HASH_FAST
		DenseHash         _map;
#endif
		KmerNb            _maxCount;
		KmerNb            _minCount;
		vector < KmerNb > _countDistribution;
		KmerNb            _nbValues;

    public:
        SimpleKmerCount ();
		void addKmer (const Kmer &kmer, bool insert=true);
		void addKmer (const Kmer &kmer, mutex &m);
		KmerNb getCount (const Kmer &kmer) const;
		bool isPresent (const KmerCode &kmerCode) const;
		bool isPresent (const Kmer &kmer) const;
		void decreaseNb (const KmerCode &kmerCode, const KmerNb nb);
		void decreaseNb (const Kmer &kmer, const KmerNb nb);
		void remove (const KmerCode &kmerCode);
		void remove (const Kmer &kmer);
		void computeCountDistribution();
		void printCountDistribution() const;
		void setMinCount(const KmerNb count);
		KmerNb getMaxCount() const;
		KmerNb getMaxCountDistribution() const;
		KmerNb getThreshold(int percent) const;
		void removeUnder(KmerNb nb);
		KmerCode getLeastFrequent();
		KmerCode getMostFrequent();
		pair <KmerCode, KmerNb> getRandom();
		bool empty() const;
		void clear();
		unsigned int getSize() const;

		friend ostream& operator<<(ostream& output, SimpleKmerCount& kc);
};

#endif



