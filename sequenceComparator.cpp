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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "sequenceComparator.hpp"


SequenceComparator::SequenceComparator(): _noMatch(false) {
	_table = new Penalty*[Globals::MAX_MERGE_SIZE+1];
	for (int i = 0; i <= Globals::MAX_MERGE_SIZE; i++) {
		_table[i] = new Penalty[Globals::MAX_MERGE_SIZE+1];
	}
}

SequenceComparator::~SequenceComparator() {
	for (int i = 0; i <= Globals::MAX_MERGE_SIZE; i++) {
		delete[] _table[i];
	}
	delete[] _table;
}

void SequenceComparator::compare (const string &s1, const string &s2) {
	_first      = s1;
	_second     = s2;
	_sizeFirst  = _first.size();
	_sizeSecond = _second.size();
	_noMatch    = false;
	initTable();
	fillTable();
	computeBackTrace();
}

void SequenceComparator::initTable () {
	int maxUnaligned = _sizeFirst - Globals::MIN_MERGE_SIZE;
	for (int i = 0; i <= Globals::MAX_MERGE_SIZE; i++) {
		_table[i][0] = (i < maxUnaligned)? i * Globals::PENALTY_SIZE: Globals::MAX_PENALTY;
		_table[0][i] = i * Globals::PENALTY_INDEL;
	}
}

void SequenceComparator::fillTable () {
	for (int i = 1; i <= _sizeFirst; i++) {
		Penalty minValue = static_cast<Penalty>(-1);
		for (int j = 1; j <= _sizeSecond; j++) {
			Penalty value1 = _table[i-1][j-1] + ((_first[i-1] == _second[j-1])? 0: Globals::PENALTY_MISMATCH);
			Penalty value2 = _table[i][j-1] + Globals::PENALTY_INDEL;
			Penalty value3 = _table[i-1][j] + Globals::PENALTY_INDEL;
			_table[i][j] = std::min<Penalty>(std::min<Penalty>(value1, value2), value3);
			minValue = std::min<Penalty>(minValue, _table[i][j]);
		}
		if (minValue >= Globals::MAX_PENALTY) {
			_noMatch = true;
			return;
		}
	}
}

void SequenceComparator::computeBackTrace () {
	if (_noMatch) {
		return;
	}
	Penalty minValue = static_cast<Penalty>(-1);
	for (int j = Globals::MIN_MERGE_SIZE + 1; j <= _sizeSecond; j++) {
		Penalty value = _table[_sizeFirst][j] + (_sizeSecond - j) * Globals::PENALTY_SIZE;
		if (value < minValue) {
			_endSecond = j;
			minValue   = value;
		}
	}
	_score = minValue;
	int i = _sizeFirst, j = _endSecond;
	int alignmentSize = 0, identity = 0;
	while ((i > _sizeFirst - Globals::MIN_MERGE_SIZE) || (j > 0)) {
		if (i == 0) {
			j--;
		}
		else if (j == 0) {
			i--;
		}
		else if (_table[i][j] == _table[i-1][j-1] + ((_first[i-1] == _second[j-1])? 0: Globals::PENALTY_MISMATCH)) {
			i--; j--;
			alignmentSize++;
			if (_first[i] == _second[j]) {
				identity++;
			}
		}
		else if (_table[i][j] == _table[i-1][j] + Globals::PENALTY_INDEL) {
			i--;
		}
		else if (_table[i][j] == _table[i][j-1] + Globals::PENALTY_INDEL) {
			j--;
		}
		else {
			cerr << "Problem in backtrace in (" << i << ", " << j << ") of table (" << _first.size() << ", " << _second.size() << "): table[" << i << "][" << j << "] = " << _table[i][j];
			if (i > 0) {
				cerr << " table[" << i-1 << "][" << j << "] = " << _table[i-1][j];
			}
			if (j > 0) {
				cerr << " table[" << i << "][" << j-1 << "] = " << _table[i][j-1];
			}
			if ((i > 0) && (j > 0)) {
				cerr << " table[" << i-1 << "][" << j-1 << "] = " << _table[i-1][j-1] << ", seq1[" << (i-1) << "] = '" << _first[i-1] << "',  seq2[" << (j-1) << "] = '" << _second[j-1] << "'";
			}
			cerr << endl;
			cerr << "Matrix is:" << endl;
			cerr << *this << endl;
		}
	}
	_startFirst = i;
	_identity   = static_cast<float>(identity) / alignmentSize;
}

Penalty SequenceComparator::getBestScore () const {
	if (_noMatch) {
		return Globals::MAX_PENALTY;
	}
	return _score;
}

float SequenceComparator::getIdentity () const {
	return _identity;
}

int SequenceComparator::getStartFirst () const {
	return _startFirst;
}

int SequenceComparator::getEndSecond () const {
	return _endSecond;
}

ostream& operator<<(ostream& output, const SequenceComparator& sq) {
	output << "\t";
	for (unsigned int i = 0; i < sq._first.size(); i++) {
		output << "\t" << sq._first[i];
	}
	output << endl;
	for (int j = 0; j <= sq._sizeSecond; j++) {
		if (j > 0)  {
			output << sq._second[j-1];
		}
		for (int i = 0; i <= sq._sizeFirst; i++) {
			output << "\t" << sq._table[i][j];
		}
		output << endl;
	}
	if (sq._noMatch) {
		output << "no match";
	}
	else {
		output << "match of " << sq._score << ", id " << sq._identity << " between " << sq._startFirst << "-" << sq._sizeFirst << " and 0-" << sq._endSecond << endl;
	}
	return output;
}

