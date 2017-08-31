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
#include "globals.hpp"
#include "kmer.hpp"


const KmerCode Kmer::UNSET;


Kmer::Kmer(): _firstCode(UNSET), _secondCode(UNSET) {}

Kmer::Kmer(const string &word): _sequence(word), _firstCode(UNSET), _secondCode(UNSET) {
	computeFirstCode();
	computeWords();
}


Kmer::Kmer(const Sequence &sequence): _sequence(sequence), _firstCode(UNSET), _secondCode(UNSET) {
	computeFirstCode();
	computeSecondCode();
}


Kmer::Kmer(const KmerCode code):  _firstCode(code), _secondCode(UNSET) {
	if (_firstCode == UNSET) {
		return;
	}
	computeSecondCode();
	computeWords();
}


Kmer::Kmer(const Kmer &k): _sequence(k._sequence), _firstCode(k._firstCode), _secondCode(k._secondCode) { }


const bool Kmer::isSet() const {
	return (_firstCode != UNSET);
}


const KmerCode Kmer::getFirstCode() const {
    return _firstCode;
}


const KmerCode Kmer::getSecondCode() const {
    return _secondCode;
}


const Sequence &Kmer::getSequence() const {
	return _sequence;
}

const string &Kmer::getFirstWord() const {
    return _sequence.getFirstWord();
}


const string &Kmer::getSecondWord() const {
    return _sequence.getSecondWord();
}


void Kmer::computeFirstCode() {
	int size = _sequence.getSize();
	_firstCode = _secondCode = 0;
    for (int i = 0; i < size; i++) {
		_firstCode  <<= Globals::NB_BITS_NUCLEOTIDES;
		_secondCode <<= Globals::NB_BITS_NUCLEOTIDES;
		_firstCode  |= Globals::getCode(_sequence.getFirstWord()[i]);
		_secondCode |= Globals::getCode(Globals::getComplement(_sequence.getFirstWord()[size - i - 1]));
	}
	if (_secondCode < _firstCode) {
		swap(_firstCode, _secondCode);
	}
}

void Kmer::computeSecondCode() {
	KmerCode firstCode = _firstCode;
	_secondCode = 0;
    for (unsigned int i = 0; i < Globals::KMER; i++) {
		_secondCode <<= Globals::NB_BITS_NUCLEOTIDES;
		_secondCode |= Globals::getComplementCode(firstCode.to_uint() & Globals::NUCLEOTIDE_MASK);
		firstCode >>= Globals::NB_BITS_NUCLEOTIDES;
	}
	if (_secondCode < _firstCode) {
		swap(_firstCode, _secondCode);
	}
}

void Kmer::computeWords() {
	KmerCode firstCode = _firstCode;
	string word  = string(Globals::KMER, 'A');
    for (unsigned int i = 0; i < Globals::KMER; i++) {
		word[Globals::KMER - i - 1] = Globals::getNucleotide(firstCode.to_uint() & Globals::NUCLEOTIDE_MASK);
		firstCode >>= Globals::NB_BITS_NUCLEOTIDES;
	}
	_sequence.setFirstWord(word);
}

const KmerCode Kmer::getCodeAfter(const int code) const {
	KmerCode newCode(code);
	newCode |= _firstCode << Globals::NB_BITS_NUCLEOTIDES;
	Kmer kmer(newCode);
	return kmer.getFirstCode();
}

const KmerCode Kmer::getCodeBefore(const int code) const {
	KmerCode newCode(code);
	Kmer kmer((newCode << (Globals::NB_BITS_NUCLEOTIDES * (Globals::KMER - 1))) | (_firstCode >> Globals::NB_BITS_NUCLEOTIDES));
	return kmer.getFirstCode();
}

const KmerCode Kmer::getCodeNeighbor(const int code) const {
	if (code > 0) {
		return getCodeAfter(code - 1);
	}
	return getCodeBefore(- code + 1);
}

const KmerCode Kmer::getCodeNeighbor(const short nucleotide, const short direction) const {
	if (direction == Globals::AFTER) {
		return getCodeAfter(nucleotide);
	}
	return getCodeBefore(nucleotide);
}

const short Kmer::compare(const Kmer &k, const short position) const {
	return (getSequence().compareNext(k.getSequence(), position));
}

ostream& operator<<(ostream& output, const Kmer& k) {
	if (! k.isSet()) {
		output << "k-mer: (empty)";
	}
	else {
		output << "k-mer " << k.getSequence() << " -- " << k._firstCode << "/" << k._secondCode;
	}
	return output;
}
