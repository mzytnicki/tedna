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

#include "globals.hpp"
#include "kmerIterator.hpp"

KmerIterator::KmerIterator(const string &s, const short size): _string(s), _size(size) {
	_currentWord = _string.substr(0, _size);
	_currentKmer = Kmer(_currentWord);
	_i           = _size;
}

void KmerIterator::getNext() {
	_i++;
	if (isOver()) {
		return;
	}
	_currentWord  = _currentWord.substr(1, _size) + _string[_i-1];
	_currentKmer  = Kmer(_currentWord);
}

bool KmerIterator::isOver() const {
	return (_i > _string.size());
}

Kmer KmerIterator::operator*() const {
	return _currentKmer;
}
