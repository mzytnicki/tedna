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
#ifndef REPEAT_HOLDER_HPP
#define REPEAT_HOLDER_HPP 1

#include <string>
#include <iostream>
#include <vector>
#include "globals.hpp"
#include "kmer.hpp"
#include "sequence.hpp"
using namespace std;

class RepeatHolder {

    private:
		Sequence         _repeat;
        KmerNb           _count;
		unsigned int     _nbKmers;
		vector<KmerCode> _kmers;

    public:
        RepeatHolder ();
        void addKmer (const Kmer &kmer, const KmerNb count);
		KmerNb getAverage () const;
		Sequence &getRepeat ();
		bool contains (const Kmer &kmer) const;
		int getNbKmers () const;
		KmerCode getKmer (const int i) const;
		void clear ();

		friend ostream& operator<<(ostream& output, const RepeatHolder& rh);
};

#endif
