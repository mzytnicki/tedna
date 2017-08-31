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

#include <algorithm>
#include "kmerSet.hpp"
#include "kmerIterator.hpp"

KmerSet::KmerSet() {}

void KmerSet::setSequence(const string &sequence) {
	_kmers.clear();
	for (unsigned int position = 0; position < sequence.size() - Globals::SHORT_KMER_SIZE + 1; position++) {
		string partF = sequence.substr(position, Globals::SHORT_KMER_SIZE);
		string partR = Globals::getReverseComplement(partF);
		string parts[] = {partF, partR};
		for (string part: parts) {
			s_kmer_t code = 0;
			for (unsigned int i = 0; i < part.size(); i++) {
				code <<= Globals::NB_BITS_NUCLEOTIDES;
				code  |= Globals::getCode(part[i]);
			}
			_kmers.insert(code);
		}
	}
}

const bool KmerSet::compare(const KmerSet &ks) const {
	for (s_kmer_t kmer1: _kmers) {
		if (ks._kmers.find(kmer1) != ks._kmers.end()) {
			return true;
		}
	}
	return false;
}

ostream& operator<<(ostream& output, const KmerSet &ks) {
	for (s_kmers_t::const_iterator it = ks._kmers.begin(); it != ks._kmers.end(); ++it) {
		output << *it << "\t";
	}
	output << endl;
	return output;
}
