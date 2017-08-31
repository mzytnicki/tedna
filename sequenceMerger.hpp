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
#ifndef SEQUENCE_MERGER_HPP
#define SEQUENCE_MERGER_HPP 1

#include <string>
#include <vector>
#include <iostream>
#include <pair>
#include "globals.hpp"
#include "sequence.hpp"
#include "repeats.hpp"
using namespace std;

class SequenceComparatorStructure {

	private:
		Penalty _penalties[2][2];
		bool    _set;

	public:
		SequenceComparatorStructure() {
			_set = false;
		}
		void setPenalty(const unsigned short position, const unsigned short sense, const Penalty penalty) {
			_set = true;
			_penalties[position][sense] = penalty;
		}
		Penalty getPenalty(const unsigned short position, const unsigned short sense) const {
			return _penalties[position][sense];
		}
		void disable() {
			_set = false;
		}
		bool isSet() const {
			return _set;
		}
};


class SequenceMerger {

    private:
		Repeats _repeats;
		int     _nbRepeats;
		SequenceComparatorStructure **_structure;

    public:
        SequenceMerger (Repeats &_repeats);
		void merge ();

	private:
		void build ();
        pair<int, int> findBestMatch () const;
		void merge(int i, int j);
		void disable(int i);
};

#endif


