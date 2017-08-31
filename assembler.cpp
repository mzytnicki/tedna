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

#include <mutex>
#include <thread>
#include "assembler.hpp"
#include "graphRepeatFinder.hpp"
#include "loopOpener.hpp"
#include "repeatMerger.hpp"
#include "inclusionRemover.hpp"
#include "scaffolder.hpp"

Assembler::Assembler(const char *fileName1, const char *fileName2, const char *outputFileName, int insertSize, int thresholdPc): _insertSize(insertSize), _thresholdPc(thresholdPc), _fileName1(fileName1), _fileName2(fileName2), _outputFileName(outputFileName) { }

void Assembler::assemble () {
	readFiles();
	computeDistributions();
	findRepeats();
	if (_repeats.empty()) {
		cout << "No repeat found. Aborting." << endl;
		exit(0);
	}
	openLoops();
	removeInclusions();
	mergeRepeats();
	scaffoldRepeats();
	removeShortRepeats();
}

void Assembler::readFiles () {
	const char *fileNames[] = {_fileName1, _fileName2};
	vector <thread> threads(Globals::NB_THREADS);
	mutex m1, m2;
	unsigned long nbReads = 0;
	for (int fileId = 0; fileId < 2; fileId++) {
		int  partId = 0;
		const char *fileName = fileNames[fileId];
		if ((Globals::NB_READS == 0) || (nbReads < Globals::NB_READS)) {
			cout << "Reading file " << (fileId+1) << "..." << endl;
			for (int threadId = 0; threadId < Globals::NB_THREADS; threadId++) {
				threads[threadId] = thread([this, fileName, &partId, &nbReads, &m1, &m2]() {
					FastqParser parser(fileName);
					unsigned long cpt = 0;
					while (! parser.isAllRead()) {
						unsigned long thisPartId;
						{
							lock_guard<mutex> lock(m1);
							thisPartId = partId;
							++partId;
							nbReads  += parser.getReadId();
						}
						if (cpt > 0) {
							cout << "\t" << nbReads << " reads read" << endl;
						}
						if ((Globals::NB_READS != 0) && (nbReads > Globals::NB_READS)) {
							cout << "\tRead enough reads." << endl;
							return;
						}
						parser.goTo(thisPartId * Globals::SIZE_THREAD, (thisPartId != 0));
						parser.endTo((thisPartId+1) * Globals::SIZE_THREAD - 1);
						for (cpt = 0; ! parser.isOver(); parser.getNextKmer(), cpt++) {
							_kmerCount.addKmer(parser.getWord(), m2);
						}
					}
				});
			}
			for (int threadId = 0; threadId < Globals::NB_THREADS; threadId++) {
				threads[threadId].join();
			}
		}
	}
}

void Assembler::computeDistributions () {
	cout << "Computing k-mer distributions..." << endl;
	_kmerCount.computeCountDistribution();
	KmerNb maxCountDistribution = _kmerCount.getMaxCountDistribution();
	cout << "\tmax count distribution: " << maxCountDistribution << " (min: " << Globals::MIN_COUNT << ")" << endl;
	cout << "\tmax count: " << _kmerCount.getMaxCount() << endl;
	if ((maxCountDistribution == Globals::MIN_COUNT) && (_thresholdPc == -1)) {
		cout << "Cannot determine maximum peak distribution. Please provide an expected repeated genome coverage." << endl;
		exit(0);
	}
	_threshold = (_thresholdPc == -1)? maxCountDistribution * Globals::NB_REPETITIONS: _kmerCount.getThreshold(_thresholdPc);
	cout << "\tthreshold " << _thresholdPc << "%: " << _threshold << endl;
	if (_threshold < Globals::MIN_COUNT) {
		cout << "\t\tincreasing threshold to " << Globals::MIN_COUNT << endl;
		_threshold = Globals::MIN_COUNT;
	}
	cout << "\tRemoving low occurrence k-mers..." << endl;
	check("Checking in the hash...");
	_kmerCount.removeUnder(_threshold);
	_kmerCount.setMinCount(_threshold);
	_kmerCount.printCountDistribution();
	check("Checking in the hash after low occurrences removal...");
}

void Assembler::findRepeats () {
	GraphRepeatFinder grf (_kmerCount, _threshold);
	grf.findRepeats();
	_repeats = grf.getRepeats();
	check("Checking in the remaining hash...");
	_kmerCount.clear();
	_repeats.check("Checking after repeat finding...");
}

void Assembler::check (const string &message) const {
	if (Globals::CHECK.empty()) {
		return;
	}
	if (Globals::CHECK.size() < Globals::KMER) {
		cout << "\t\tWarning! Cannot check given sequence for its size is less than k-mer size." << endl;
		return;
	}
	cout << "\t\t" << message << endl;
	for (unsigned int position = 0; position < Globals::CHECK.size() - Globals::KMER + 1; position++) {
		string part  = Globals::CHECK.substr(position, Globals::KMER);
		KmerNb count = _kmerCount.getCount(Kmer(part));
		if (count > 0) {
			cout << "\t\t\tGot '" << part << "' @ " << position << ", " << count << " times" << endl;
		}
	}
}

void Assembler::openLoops () {
	//cout << "Opening loops with\n" << _repeats << endl;
	LoopOpener loopOpener (_fileName1, _repeats);
	loopOpener.openLoops();
	_repeats = loopOpener.getRepeats();
	_repeats.check("Checking after loop opening...");
}

void Assembler::removeInclusions () {
	//cout << "Removing inclusions with\n" << _repeats << endl;
	if (_repeats.getNbRepeats() > 1) {
		InclusionRemover ir (_repeats);
		ir.removeInclusions();
		_repeats = ir.getRepeats();
		_repeats.check("Checking after inclusion removal...");
	}
}

void Assembler::mergeRepeats () {
	//cout << "Merging repeats with\n" << _repeats << endl;
	if (_repeats.getNbRepeats() > 1) {
		RepeatMerger rm (_repeats, _threshold);
		rm.mergeRepeats();
		_repeats = rm.getRepeats();
		_repeats.check("Checking after repeat merging...");
	}
}

void Assembler::scaffoldRepeats () {
	//cout << "Scaffolding repeats with\n" << _repeats << endl;
	if (_repeats.getNbRepeats() > 1) {
		Scaffolder s (_repeats, _fileName1, _fileName2, _insertSize);
		s.scaffold();
		_repeats = s.getRepeats();
		_repeats.check("Checking after scaffolding...");
	}
}

void Assembler::removeShortRepeats () {
	for (unsigned int i = 0; i < _repeats.getNbRepeats(); i++) {
		if ((_repeats[i].getSize() < Globals::MIN_TE_SIZE) || (_repeats[i].getSize() > Globals::MAX_TE_SIZE)) {
			_repeats.removeRepeat(i);
		}
	}
}

void Assembler::dump () {
	ofstream handle;
	handle.open(_outputFileName);
	_repeats.printFasta(handle);
	handle.close();
	cout << "Found " << _repeats.getNbRepeats() << " putative TEs." << endl;
}
