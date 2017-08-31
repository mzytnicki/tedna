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
#include <thread>
#include "loopOpener.hpp"
#include "kmerIterator.hpp"
#include "fastqParser.hpp"

LoopOpener::LoopOpener(const char *fileName, const Repeats &repeats): _fileName(fileName), _repeats(repeats), _nbLoops(0) { }

void LoopOpener::openLoops() {
	for (unsigned int i = 0; i < _repeats.getNbRepeats(); i++) {
		if (_repeats[i].isLoop()) {
			addLoop(_repeats[i].getRepeat());
			_nbLoops++;
		}
	}
	if (empty()) {
		return;
	}
	readReads();
	for (unsigned int i = 0; i < _repeats.getNbRepeats(); i++) {
		if (_repeats[i].isLoop()) {
			string sequence        = _repeats[i].getRepeat().getFirstWord();
			pair <int, int> region = getOverRepresented(sequence);
			//cout << "Region is " << region.first << ", " << region.second << endl;
			if (region.first != -1) {
				string newSequence = sequence.substr(region.first) + ((region.first <= region.second)? "": sequence.substr(Globals::KMER-1)) + sequence.substr(Globals::KMER-1, region.second+1);
				//cout << "size from " << sequence.size() << " to " << newSequence.size() << " with count " << _pathes[loopId].getCount() << endl;
				Sequence s(newSequence);
				_repeats.changeRepeat(i, s, _repeats[i].getCount());
				_repeats.setLoop(i, false);
			}
		}
	}
}

const Repeats &LoopOpener::getRepeats() {
	return _repeats;
}

void LoopOpener::reset() {
	_count.clear();
}

void LoopOpener::addLoop(const Sequence &loop) {
	KmerCode currentCode, otherCode;
	int i = 0;
	for (KmerIterator it = KmerIterator(loop.getFirstWord()); ! it.isOver(); it.getNext(), i++) {
		currentCode = (*it).getFirstCode();
		_count[currentCode] = 0;
	}
}

void LoopOpener::readReads() {
	cout << "Starting loop opener, " << _nbLoops << " loop(s) found..." << endl;
	vector <thread> threads(Globals::NB_THREADS);
	mutex m1, m2, m3;
	int   partId = 0;
	unsigned long nbReads = 0;
	for (int threadId = 0; threadId < Globals::NB_THREADS; threadId++) {
		threads[threadId] = thread([this, &partId, &nbReads, &m1, &m2, &m3]() {
			FastqParser parser(_fileName);
			Kmer currentKmer;
			KmerCode currentCode;
			while (! parser.isAllRead()) {
				long unsigned thisPartId;
				if (nbReads > 0) {
					cout << "\t" << nbReads << " reads read." << endl;
				}
				if ((Globals::NB_READS != 0) && (nbReads > Globals::NB_READS)) {
					cout << "\t" << nbReads << " reads read." << endl;
					return;
				}
				{
					lock_guard<mutex> lock(m1);
					thisPartId = partId;
					nbReads   += parser.getReadId();
					++partId;
				}
				parser.goTo(thisPartId * Globals::SIZE_THREAD, (thisPartId != 0));
				parser.endTo((thisPartId+1) * Globals::SIZE_THREAD - 1);
				for (unsigned long i = 0; !parser.isOver(); i++, parser.getNextKmer()) {
					currentKmer = Kmer(parser.getSequence());
					currentCode = currentKmer.getFirstCode();
					if (_count.find(currentCode) != _count.end()) {
						lock_guard<mutex> lock(m2);
						_count[currentCode]++;
					}
				}
			}
		});
	}
	for (int threadId = 0; threadId < Globals::NB_THREADS; threadId++) {
		threads[threadId].join();
	}
}


pair <int, int> LoopOpener::getOverRepresented(const string &loop) {
	vector < pair <int, int> > regions, mergedRegions;
	KmerNb sum = 0, avg, count;
	int i = 0, start = 0, biggest = 0, loopSize, ltrSize;
	pair <int, int> bestRegion;
	bool inRegion = false;
	for (map<KmerCode, KmerNb>::iterator it = _count.begin(); it != _count.end(); ++it) {
		sum += it->second;
	}
	avg = sum / _count.size();
	//cout << "threshold: " << avg * 1.5 << endl;
	for (KmerIterator it = KmerIterator(loop); ! it.isOver(); it.getNext(), ++i) {
		count = _count[(*it).getFirstCode()];
		//cout << "position: " << i << ": " << count << endl;
		if ((count > 1.5 * avg) && (! inRegion)) {
			start    = i;
			inRegion = true;
		}
		else if ((count <= 1.5 * avg) && (inRegion)) {
			regions.push_back(pair<int, int>(start, i-1));
			inRegion = false;
		}
	}
	loopSize = i-1;
	regions.push_back(pair<int, int>(start, loopSize));
	start = regions.front().first;
	//cout << "region: " << regions.front().first << ", " << regions.front().second << endl;
	for (unsigned int i = 1; i < regions.size(); i++) {
		//cout << "region: " << regions[i].first << ", " << regions[i].second << endl;
		if (regions[i].first - regions[i-1].second > 10) {
			mergedRegions.push_back(pair<int, int>(start, regions[i-1].second));
			start = regions[i].first;
		}
	}
	mergedRegions.push_back(pair<int, int>(start, regions.back().second));
	for (unsigned int i = 0; i < mergedRegions.size(); i++) {
		//cout << "merged region: " << mergedRegions[i].first << ", " << mergedRegions[i].second << endl;
		if (mergedRegions[i].second - mergedRegions[i].first > biggest) {
			bestRegion = mergedRegions[i];
			biggest = bestRegion.second - bestRegion.first;
		}
	}
	ltrSize = bestRegion.second - bestRegion.first + 1;
	if ((mergedRegions.front().first == 0) && (mergedRegions.back().second == loopSize) && (mergedRegions.front().second + loopSize - mergedRegions.back().first > biggest)) {
		bestRegion.first  = mergedRegions.back().first;
		bestRegion.second = mergedRegions.front().second;
		ltrSize           = loopSize - bestRegion.first + bestRegion.second + 1;
	}
	//cout << "best region: " << bestRegion.first << ", " << bestRegion.second << endl;
	if ((ltrSize < Globals::MIN_LTR_SIZE) || (ltrSize > Globals::MAX_LTR_SIZE)) {
		return make_pair(-1, -1);
	}
	return bestRegion;
}

bool LoopOpener::empty() const {
	return (_count.empty());
}

ostream& operator<<(ostream& output, const LoopOpener& lo) {
	for (map < KmerCode, KmerNb >::const_iterator it = lo._count.begin(); it != lo._count.end(); it++) {
		output << Kmer(it->first) << " (" << it->second << ")   ";
	}
	return output;
}

