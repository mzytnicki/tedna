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
#ifndef KMER_PARSER_HPP
#define KMER_PARSER_HPP 1

#include <fstream>
#include "globals.hpp"
#include "hashes.hpp"
#include "countDistribution.hpp"
using namespace std;

class KmerParser {

    protected:
	    ifstream _file;
	    KmerNb   _maxCount;
	    KmerNb   _nbValues;

    public:
	    KmerParser (const char *fileName);
	    void getStats ();
	    KmerNb getMaxCount () const;
	    KmerNb getNbValues () const;
	    void fill (CountDistribution &cd);
	    void fill (hash_t &map, KmerNb minCount);
	    //void fillKmers (BBhashStr &bb, KmerNb minCount);
	    //void fillCounts (BBhashStr &bb, KmerNb minCount);
};

#endif
