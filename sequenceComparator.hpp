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
#ifndef SEQUENCE_COMPARATOR_HPP
#define SEQUENCE_COMPARATOR_HPP 1

#include <string>
#include <iostream>
#include "globals.hpp"
using namespace std;


class SequenceComparator {

    private:
		//Penalty _table[Globals::MAX_MERGE_SIZE+1][Globals::MAX_MERGE_SIZE+1];
		Penalty **_table;
		bool    _noMatch;
		string  _first, _second;
		int     _sizeFirst, _sizeSecond;
		int     _startFirst, _endSecond;
		Penalty _score;
		float   _identity;

    public:
        SequenceComparator ();
        ~SequenceComparator ();
		void compare (const string &s1, const string &s2);
		Penalty getBestScore () const;
		float getIdentity () const;
		int getStartFirst () const;
		int getEndSecond () const;

		friend ostream& operator<<(ostream& output, const SequenceComparator& sq);

	private:
		void initTable ();
		void fillTable ();
		void computeBackTrace ();

};

#endif
