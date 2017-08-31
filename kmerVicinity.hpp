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
#ifndef KMER_VICINITY_HPP
#define KMER_VICINITY_HPP 1

#include <string>
#include <iostream>
#include "globals.hpp"
#include "kmer.hpp"
using namespace std;

class KmerVicinity {

    private:
		unsigned int _count;
		unsigned int _countBefore[Globals::NB_NUCLEOTIDES];
		unsigned int _countAfter[Globals::NB_NUCLEOTIDES];

    public:
        KmerVicinity ();
		void addCount ();
		void addCountBefore (const int i);
		void addCountAfter (const int i);
        const unsigned int getCount  () const;
        const unsigned int getCountBefore (const int i) const;
        const unsigned int getCountAfter (const int i) const;
        const int getMaxCountBefore () const;
        const int getMaxCountAfter () const;

		friend ostream& operator<<(ostream& output, const KmerVicinity& kv);
};

#endif
