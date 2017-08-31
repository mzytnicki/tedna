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
#ifndef KMERSET_HPP
#define KMERSET_HPP 1

#include <string>
#include <set>
#include <iostream>
using namespace std;

typedef unsigned int s_kmer_t;
typedef set < s_kmer_t > s_kmers_t;

class KmerSet {

    private:
		s_kmers_t _kmers;

    public:
        KmerSet ();

		void setSequence (const string &sequence);
		const bool compare (const KmerSet &ks) const;

		friend ostream& operator<<(ostream& output, const KmerSet& ks);
};

#endif
