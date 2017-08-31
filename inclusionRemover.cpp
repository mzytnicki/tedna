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

#include <thread>
#include <mutex>
#include "globals.hpp"
#include "inclusionRemover.hpp"

InclusionRemover::InclusionRemover(const Repeats &repeats): _kmerSets(repeats.getNbRepeats()), _repeats(repeats) { }

void InclusionRemover::removeInclusions () {
	cout << "Removing duplicates (" << _repeats.getNbRepeats() << " elements)..." << endl;
	_repeats.sort();
	storeKmers();
	unsigned int i = 0, j = 1, cpt = 0, nbInclusions = 0;
	vector <thread> threads(Globals::NB_THREADS);
	mutex m1, m2;
	for (int threadId = 0; threadId < Globals::NB_THREADS; threadId++) {
		threads[threadId] = thread([this, &i, &j, &cpt, &nbInclusions, &m1, &m2]() {
			unsigned int thisI, thisJ, thisCpt, nbRepeats = _repeats.getNbRepeats();
			while (true) {
				{
					lock_guard<mutex> lock(m1);
					thisI = i;
					thisJ = j;
					do {
						j++;
						if (j >= nbRepeats) {
							i++;
							j = i+1;
						}
					}
					while ((j < nbRepeats) && ((_repeats.isRemoved(i)) || (_repeats.isRemoved(j))));
					cpt++;
					thisCpt = cpt;
				}
				if (thisJ >= nbRepeats) {
					return;
				}
				if ((! _repeats.isRemoved(thisI)) && (! _repeats.isRemoved(thisJ))) {
					if (checkInclusion(thisI, thisJ)) {
						lock_guard<mutex> lock(m2);
						nbInclusions++;
						_repeats.removeRepeat(thisJ);
					}
				}
				if (thisCpt % 100000 == 0) {
					cout << "\t" << thisCpt << " duplications evaluated." << endl;
				}
			}
		});
	}
	for (int threadId = 0; threadId < Globals::NB_THREADS; threadId++) {
		threads[threadId].join();
	}
	cout << "\t" << cpt << " duplications evaluated, " << nbInclusions << " found." << endl;
	_repeats.sort();
	//cout << *this << endl;
}

void InclusionRemover::storeKmers () {
	for (unsigned int i = 0; i < _repeats.getNbRepeats(); i++) {
		if (! _repeats.isRemoved(i)) {
			_kmerSets[i].setSequence(_repeats[i].getRepeat().getFirstWord());
			//cout << "\tStoring kmers[" << i << "]:" << _kmerSets[i] << endl;
		}
	}
}

bool InclusionRemover::checkInclusion (const unsigned int i, const unsigned int j) const {
	//cout << "\tInspecting " << i << " vs " << j << endl;
	if (! _kmerSets[i].compare(_kmerSets[j])) {
		//cout << "\t\tNo kmer in common." << endl;
		return false;
	}
	//cout << "\t\tSome kmers in common." << endl;
	const string &firstString = _repeats.getRepeat(i).getRepeat().getFirstWord();
	for (short int direction = 0; direction < Globals::DIRECTIONS; direction++) {
		const string &secondString = _repeats.getRepeat(j).getRepeat().getWord(direction);
		bool   included            = compareStrings(secondString, firstString);
		if (included) {
			//cout << secondString << " seems included into " << firstString << endl;
			return true;
		}
	}
	return false;
}

bool InclusionRemover::compareStrings(const string &firstString, const string &secondString) const {
	//cout << "Entering compareStrings" << endl;
	unsigned int firstSize = firstString.size(), secondSize = secondString.size();
	vector < int > line(firstSize+1);
	int current, tmp;
	int maxPenalty = firstSize * Globals::MAX_IDENTITY;
	//cout << "\t\tfirst string: " << firstString << ", second string: " << secondString << ", max penalty: " << maxPenalty << endl;
	for (unsigned int i = 1; i <= firstSize; i++) {
		line[i] = i;
	}
	for (unsigned int i = 1; i <= secondSize; i++) {
		current = 0;
		for (unsigned int j = 1; j <= firstSize; j++) {
			tmp = line[j];
			try {
				line[j] = min<int>(min<int>(current + ((secondString[i-1]==firstString[j-1])? 0: 1), line[j-1]+1), line[j]+1);
			}
			catch (...) {
				return false;
			}
			current = tmp;
			//cout << table[i][j] << " ";
		}
		if (line[firstSize] < maxPenalty) {
			//cout << "\t\tfound good hit at (" << i << ", " << secondString.size() << ") with penalty " << table[i][secondString.size()] << endl;
			//cout << "Leaving it 1" << endl;
			return true;
		}
		//cout << endl;
	}
	//cout << "Leaving it 2" << endl;
	return false;
}

const Repeats &InclusionRemover::getRepeats () {
	return _repeats;
}

ostream& operator<<(ostream& output, const InclusionRemover& rm) {
	output << "Sequences:\n" << rm._repeats;
	return output;
}
