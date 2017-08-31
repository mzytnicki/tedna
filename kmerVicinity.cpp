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

#include "kmerVicinity.hpp"

KmerVicinity::KmerVicinity(): _count(0) {
	for (int i = 0; i < Globals::NB_NUCLEOTIDES; i++) {
		_countBefore[i] = _countAfter[i] = 0;
	}
}


void KmerVicinity::addCount() {
	if (_count != static_cast<unsigned int>(-1)) {
		_count++;
	}
}


void KmerVicinity::addCountBefore(const int i) {
	if (_countBefore[i] != static_cast<unsigned int>(-1)) {
		_countBefore[i]++;
	}
}


void KmerVicinity::addCountAfter(const int i) {
	if (_countAfter[i] != static_cast<unsigned int>(-1)) {
		_countAfter[i]++;
	}
}


const unsigned int KmerVicinity::getCount() const {
    return _count;
}


const unsigned int KmerVicinity::getCountBefore(const int i) const {
    return _countBefore[i];
}


const unsigned int KmerVicinity::getCountAfter(const int i) const {
    return _countAfter[i];
}


const int KmerVicinity::getMaxCountBefore() const {
	unsigned int value = 0;
	int index = 0;
	for (int i = 0; i < Globals::NB_NUCLEOTIDES; i++) {
		if (value > _countBefore[i]) {
			value = _countBefore[i];
			index = i;
		}
	}
	return index;
}


const int KmerVicinity::getMaxCountAfter() const {
	unsigned int value = 0;
	int index = 0;
	for (int i = 0; i < Globals::NB_NUCLEOTIDES; i++) {
		if (value > _countAfter[i]) {
			value = _countAfter[i];
			index = i;
		}
	}
	return index;
}


ostream& operator<<(ostream& output, const KmerVicinity& k) {
	for (int i = 0; i < Globals::NB_NUCLEOTIDES; i++) {
		output << " " << Globals::getNucleotide(i) << ": " << k.getCountBefore(i) << " -";
	}
	output << "- " << k.getCount() << " -";
	for (int i = 0; i < Globals::NB_NUCLEOTIDES; i++) {
		output << "- " << Globals::getNucleotide(i) << ": " << k.getCountAfter(i) << " ";
	}
	return output;
}
