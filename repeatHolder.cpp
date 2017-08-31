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

#include "repeatHolder.hpp"

RepeatHolder::RepeatHolder(): _count(0), _nbKmers(0) { }


void RepeatHolder::addKmer(const Kmer &kmer, const KmerNb count) {
	_count   += count;
	_nbKmers ++;
	_kmers.push_back(kmer.getFirstCode());
	if (_repeat.empty()) {
		_repeat = kmer.getFirstWord();
		return;
	}
	int comparison = _repeat.compareNext(kmer.getSequence());
	_repeat.addCode(comparison);
}


KmerNb RepeatHolder::getAverage () const {
	return static_cast<KmerNb>(_count / _nbKmers);
}


Sequence &RepeatHolder::getRepeat() {
	return _repeat;
}


bool RepeatHolder::contains(const Kmer &kmer) const {
	KmerCode code = kmer.getFirstCode();
	for (unsigned int i = 0; i < _kmers.size(); i++) {
		if (code == _kmers[i]) {
			return true;
		}
	}
	return false;
}


int RepeatHolder::getNbKmers () const {
	return _kmers.size();
}


KmerCode RepeatHolder::getKmer (const int i) const {
	return _kmers[i];
}


void RepeatHolder::clear () {
	_repeat.clear();
	_count   = 0;
	_nbKmers = 0;
	_kmers.clear();
}


ostream& operator<<(ostream& output, const RepeatHolder& rh) {
	output << rh._repeat << " (" << rh.getAverage() << ")";
	return output;
}
