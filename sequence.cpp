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

#include <cmath>
#include <sstream>
#include "globals.hpp"
#include "sequence.hpp"


Sequence::Sequence() {}

Sequence::Sequence(const string &word): _firstWord(word) {
	computeSecondWord();
}


Sequence::Sequence(const Sequence &s): _firstWord(s._firstWord), _secondWord(s._secondWord) { }


const bool Sequence::empty() const {
	return _firstWord.empty();
}


const unsigned int Sequence::getSize() const {
	return _firstWord.size();
}

const unsigned int Sequence::getUnambiguousSize() const {
	unsigned int size = 0;
	for (char c: _firstWord) {
		switch (c) {
			case 'A':
			case 'C':
			case 'G':
			case 'T':
			case 'a':
			case 'c':
			case 'g':
			case 't':
				size++;
		}
	}
	return size;
}

const string &Sequence::getFirstWord() const {
    return _firstWord;
}


const string &Sequence::getSecondWord() const {
    return _secondWord;
}


const string &Sequence::getWord(const short i) const {
	return (i == Globals::DIRECT)? _firstWord: _secondWord;
}


void Sequence::setFirstWord(const string &word) {
	_firstWord = word;
	computeSecondWord();
}


bool Sequence::addFront(const string &nucleotides) {
	_firstWord.insert(0, nucleotides);
	return (computeSecondWord());
}


bool Sequence::addBack(const string &nucleotides) {
	_firstWord.append(nucleotides);
	return (computeSecondWord());
}


bool Sequence::addCode(const int code) {
	if (code < 0) {
		addFront(string(1, Globals::getNucleotide(-code-1)));
		return updateWords();
	}
	if (code > 0) {
		addBack(string(1, Globals::getNucleotide(code-1)));
		return updateWords();
	}
	cerr << "Cannot use code 0 for " << *this << endl;
	return false;
}


void Sequence::clear() {
	if (! empty()) {
		_firstWord.clear();
		_secondWord.clear();
	}
}


bool Sequence::updateWords() {
	if (_secondWord < _firstWord) {
		swap(_firstWord, _secondWord);
		return true;
	}
	return false;
}


bool Sequence::computeSecondWord() {
	_secondWord = Globals::getReverseComplement(_firstWord);
	return updateWords();
}


const short Sequence::compareNext(const Sequence &k, const short position) const {
	unsigned int thisSize = getSize(), thatSize = k.getSize();
	int last = min<unsigned int>(thisSize, thatSize) - 1;
	if ((position == Globals::AFTER) || (position == Globals::POSITIONS)) {
		if (k.getFirstWord().substr(0, last) == getFirstWord().substr(thisSize-last, last)) {
			if (position == Globals::AFTER) {
				return Globals::DIRECT;
			}
			return Globals::getCode(k.getFirstWord()[last]) + 1;
		}
		if (k.getSecondWord().substr(0, last) == getFirstWord().substr(thisSize-last, last)) {
			if (position == Globals::AFTER) {
				return Globals::REVERSE;
			}
			return Globals::getCode(k.getSecondWord()[last]) + 1;
		}
	}
	else if ((position == Globals::BEFORE) || (position == Globals::POSITIONS)) {
		if (k.getFirstWord().substr(thatSize-last, last) == getFirstWord().substr(0, last)) {
			if (position == Globals::BEFORE) {
				return Globals::DIRECT;
			}
			return -(Globals::getCode(k.getFirstWord()[0]) + 1);
		}
		if (k.getSecondWord().substr(thatSize-last, last) == getFirstWord().substr(0, last)) {
			if (position == Globals::BEFORE) {
				return Globals::REVERSE;
			}
			return -(Globals::getCode(k.getSecondWord()[0]) + 1);
		}
	}
	return 0;
}


tuple<short, short, bool> Sequence::merge(const Sequence &k, short direction, int size) {
	short AFTER = Globals::AFTER, BEFORE = Globals::BEFORE, DIRECT = Globals::DIRECT, REVERSE = Globals::REVERSE;
	tuple<short, short, bool> code = make_tuple(AFTER, DIRECT, false);
	if (empty()) {
		*this = k;
		return code;
	}
	if (k.empty()) {
		return code;
	}
	unsigned int thisSize = getSize(), thatSize = k.getSize();
	int maxSize  = min<unsigned int>(thisSize, thatSize) - 1;
	int minSize  = 5;
	if (size != -1) {
		minSize = size;
		maxSize = size;
	}
	//cout << "\tMerging " << *this << " with " << k << " and direction " << direction;
	for (int size = maxSize; size >= minSize; size--) {
		if (direction != Globals::BEFORE) {
			if (k.getFirstWord().substr(0, size) == getFirstWord().substr(thisSize-size, size)) {
				//cout << " is " << *this << " (size: " << thisSize << ", " << thatSize << ", " << size << ", case 1)" << endl;
				bool reverse = addBack(k.getFirstWord().substr(size, thatSize-size));
				return make_tuple(AFTER, DIRECT, reverse);
			}
			if (k.getSecondWord().substr(0, size) == getFirstWord().substr(thisSize-size, size)) {
				//cout << " is " << *this << " (size: " << thisSize << ", " << thatSize << ", " << size << ", case 2)" << endl;
				bool reverse = addBack(k.getSecondWord().substr(size, thatSize-size));
				return make_tuple(AFTER, REVERSE, reverse);
			}
		}
		if (direction != Globals::AFTER) {
			if (k.getFirstWord().substr(thatSize-size, size) == getFirstWord().substr(0, size)) {
				//cout << " is " << *this << " (size: " << thisSize << ", " << thatSize << ", " << size << ", case 3)" << endl;
				bool reverse = addFront(k.getFirstWord().substr(0, thatSize-size));
				return make_tuple(BEFORE, DIRECT, reverse);
			}
			if (k.getSecondWord().substr(thatSize-size, size) == getFirstWord().substr(0, size)) {
				//cout << " is " << *this << " (size: " << thisSize << ", " << thatSize << ", " << size << ", case 4)" << endl;
				bool reverse = addFront(k.getSecondWord().substr(0, thatSize-size));
				return make_tuple(BEFORE, REVERSE, reverse);
			}
		}
	}
	cerr << "Cannot merge " << *this << " with " << k << " with direction " << direction << endl;
	exit(1);
	return make_tuple(-1, -1, false);
}

bool Sequence::isAmbiguous () const {
	for (char c: _firstWord) {
		switch(c) {
			case 'A':
			case 'a':
			case 'C':
			case 'c':
			case 'G':
			case 'g':
			case 'T':
			case 't':
			case 'U':
			case 'u':
				break;
			default:
				return true;
		}
	}
	return false;
}

bool Sequence::isLowComplexity () const {
	bool  nucleotides[Globals::NB_NUCLEOTIDES];
	short nbNucleotides = 0;
	for (int i = 0; i < Globals::NB_NUCLEOTIDES; i++) {
		nucleotides[i] = false;
	}
	for (unsigned int i = 0; i < _firstWord.length(); i++) {
		short code = Globals::getCode(_firstWord[i]);
		if (! nucleotides[code]) {
			nbNucleotides++;
			if (nbNucleotides >= 3) {
				return false;
			}
			nucleotides[code] = true;
		}
	}
	return true;
}

string Sequence::printFasta(string title) const {
	stringstream fastaString;
	string word     = _firstWord;
	int    lineSize = 60;
	fastaString << ">" << title << " (" << _firstWord.size() << ")";
	for (unsigned int i = 0; i <= _firstWord.length()/lineSize; i++) {
		fastaString << "\n" << _firstWord.substr(i*lineSize, lineSize);
	}
	fastaString << "\n";
	return fastaString.str();
}

bool operator==(const Sequence &s1, const Sequence &s2) {
	return (s1._firstWord == s2._firstWord);
}

ostream& operator<<(ostream& output, const Sequence& s) {
	output << s._firstWord << "/" << s._secondWord;
	return output;
}
