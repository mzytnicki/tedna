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
#ifndef REDUCED_KMER_COUNT_HPP
#define REDUCED_KMER_COUNT_HPP 1

#include <cmath>
#include <string>
#include <iostream>
#include "kmer.hpp"
using namespace std;

class KmerVicinity {

};

class ReducedKmerCount {

    private:
		int             _size;
		unsigned int    _countSize;
        unsigned int    _number;
		unsigned int*   _mappingTable;
		unsigned int*   _countDistribution;

    public:
        KmerCount (const int size);
        ~KmerCount ();
        void addKmer (const Kmer &kmer);
        int getCount (const Kmer &kmer);
		void computeCountDistribution ();
		void printCountDistribution () const;
		unsigned int getMaxCountDistribution () const;

		friend ostream& operator<<(ostream& output, const KmerCount& kc);
};

#endif

