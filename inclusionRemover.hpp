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
#ifndef INCLUSION_REMOVER_HPP
#define INCLUSION_REMOVER_HPP 1

#include <vector>
#include <iostream>
#include "globals.hpp"
#include "repeats.hpp"
#include "kmerSet.hpp"
using namespace std;



class InclusionRemover {

    private:
		vector <KmerSet> _kmerSets;
		Repeats  _repeats;
		int      _i, _j;
		short    _direction, _position;
		KmerNb   _minCount;

    public:
        InclusionRemover (const Repeats &repeats);
        void removeInclusions ();
		const Repeats &getRepeats();

	private:
		void storeKmers ();
		bool checkInclusion (const unsigned int i, const unsigned int j) const;
		bool compareStrings (const string &firstString, const string &secondString) const;

		friend ostream& operator<<(ostream& output, const InclusionRemover& ir);
};

#endif
