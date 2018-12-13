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

#include <map>
#include <tuple>
#include <limits>
#include <algorithm>
#include "globals.hpp"
#include "scaffolder.hpp"
#include "graphTrimmer.hpp"


Scaffolder::Scaffolder (Repeats &r, const char *fileName1, const char *fileName2, unsigned int insertSize): _inputRepeats(r), _fileName1(fileName1), _fileName2(fileName2), _insertSize(insertSize), _maxEvidencesPerNode(1) { }

void Scaffolder::scaffold () {
	cout << "Starting scaffolding (" << _inputRepeats.getNbRepeats() << " fragments)..." << endl;
	//cout << _inputRepeats << endl;
	storeKmers();
	cout << "\tRepeats read..." << endl;
	//buildStructure();
	//cout << "\tStructure built..." << endl;
	fillStructure();
	cout << "\tStructure filled..." << endl;
	//removeWeakLinks();
	//cout << "\tWeak links removed..." << endl;
	//computeMode();
	//cout << "\tCounts checked..." << endl;
	//checkDistances();
	//cout << "\tDistances checked..." << endl;
	if (Globals::SCAFFOLD_MAX_EV != 0) {
		resetMinEvidences();
		if (_maxEvidencesPerNode > 1) {
			cout << "\tMin # of evidences adjusted to " << _maxEvidencesPerNode << endl;
		}
	}
	buildGraphs();
}

void Scaffolder::storeKmers () {
	for (unsigned int repeatId = 0; repeatId < _inputRepeats.getNbRepeats(); repeatId++) {
		string repeat = _inputRepeats.getRepeat(repeatId).getRepeat().getFirstWord();
		for (unsigned int position = 0; position < repeat.size() - Globals::KMER + 1; position++) {
			string partF = repeat.substr(position, Globals::KMER);
			if (! Sequence(partF).isAmbiguous()) {
				string partR = Globals::getReverseComplement(partF);
				string parts[] = {partF, partR};
				for (short direction = 0; direction < Globals::DIRECTIONS; direction++) {
					_kmers[parts[direction]].push_back(tuple<int, int, short>(repeatId, position, direction));
					//cout << "repeat kmers: # repeat " << repeatId << ", position: " << position << ", direction: " << direction << " --> " << parts[direction] << endl;
				}
			}
		}
	}
	/*
	for (pair <const string, vector < tuple <int, int, short> > > i: _kmers) {
		cout << i.first << " ->" << endl;
		for (tuple <int, int, short > t: i.second) {
			cout << "\trepeatId: " << get<0>(t) << ", position: " << get<1>(t) << ", direction: " << get<2>(t) << endl;
		}
	}
	*/
}

/*
void Scaffolder::buildStructure () {
	_distances = vector < vector < array < array < vector <int>, Globals::DIRECTIONS >, Globals::POSITIONS > > >(_inputRepeats.getNbRepeats());
	for (unsigned int i = 1; i < _inputRepeats.getNbRepeats(); i++) {
		_distances[i] = vector < array < array < vector <int>, Globals::DIRECTIONS >, Globals::POSITIONS > >(i);
	}
	_modes = vector < vector < array < array < int, Globals::DIRECTIONS >, Globals::POSITIONS > > >(_inputRepeats.getNbRepeats());
	for (unsigned int i = 1; i < _inputRepeats.getNbRepeats(); i++) {
		_modes[i] = vector < array < array < int, Globals::DIRECTIONS >, Globals::POSITIONS > >(i);
		for (unsigned int j = 0; j < i; j++) {
			for (int k = 0; k < Globals::POSITIONS; k++) {
				for (int l = 0; l < Globals::DIRECTIONS; l++) {
					_modes[i][j][k][l] = numeric_limits<int>::min();
					_distances[i][j][k][l].clear();
				}
			}
		}
	}
}
*/

void Scaffolder::fillStructure () {
	vector <thread> threads(Globals::NB_THREADS);
	mutex m1, m2;
	int   partId = 0;
	unsigned long nbReads = 0;
	for (int threadId = 0; threadId < Globals::NB_THREADS; threadId++) {
		threads[threadId] = thread([this, &partId, &nbReads, &m1, &m2]() {
			FastxParser *parser1, *parser2;
			if (Globals::FASTA_INPUT) {
				parser1 = new FastaParser(_fileName1);
				parser2 = new FastaParser(_fileName2);
			}
			else {
				parser1 = new FastqParser(_fileName1);
				parser2 = new FastqParser(_fileName2);
			}
			parser1->reset();
			while (! parser1->isAllRead()) {
				unsigned long thisPartId;
				if (nbReads > 0) {
					cout << "\t" << nbReads << " pairs read." << endl;
				}
				if ((Globals::NB_READS != 0) && (nbReads > Globals::NB_READS)) {
					cout << "\tRead enough pairs." << endl;
					return;
				}
				{
					lock_guard<mutex> lock(m1);
					thisPartId = partId;
					nbReads   += parser1->getReadId();
					++partId;
				}
				parser1->goTo(thisPartId * Globals::SIZE_THREAD, (partId != 0));
				parser1->endTo((thisPartId+1) * Globals::SIZE_THREAD - 1);
				parser2->goTo(thisPartId * Globals::SIZE_THREAD, (partId != 0));
				parser2->endTo((thisPartId+1) * Globals::SIZE_THREAD - 1);
				fillStructure(parser1, parser2, m2);
			}
			delete parser1;
			delete parser2;
		});
	}
	for (int threadId = 0; threadId < Globals::NB_THREADS; threadId++) {
		threads[threadId].join();
	}
	//cout << "done" << endl;
}

void Scaffolder::fillStructure (FastxParser *parser1, FastxParser *parser2, mutex &m) {
	//cout << _inputRepeats << endl;
	for (; ! parser1->isOver(); parser1->getNextLine(), parser2->getNextLine()) {
		//cout << "lines: " << parser1.getLine() << "\t" << parser2.getLine() << endl;
		//string lines[] = {parser1.getLine(), parser2.getLine()};
		string lines[] = {parser1->getLine(), Globals::getReverseComplement(parser2->getLine())};
		if ((lines[0].size() >= Globals::KMER) && (lines[1].size() >= Globals::KMER)) {
			vector < tuple < int, int, int, short > > positions2;
			//cout << "Positions2" << endl;
			for (unsigned int readPos2 = 0; readPos2 < lines[1].size() - Globals::KMER + 1; readPos2++) {
				string part  = lines[1].substr(readPos2, Globals::KMER);
				auto   match = _kmers.find(part);
				if (match != _kmers.end()) {
					for (tuple <int, int, short> t: match->second) {
						positions2.push_back(make_tuple(readPos2, get<0>(t), get<1>(t), get<2>(t)));
						//cout << "\tword: " << lines[1].substr(readPos2, Globals::KMER) << ", pos: " << readPos2 << ", repeat #: " << get<0>(t) << ", repeat pos: " << get<1>(t) << ", repeat strand: " << get<2>(t) << endl;
					}
				}
			}
			//cout << "Positions1" << endl;
			for (unsigned int readPos1 = 0; readPos1 < lines[0].size() - Globals::KMER + 1; readPos1++) {
				string part  = lines[0].substr(readPos1, Globals::KMER);
				auto   match = _kmers.find(part);
				if (match != _kmers.end()) {
					for (tuple <int, int, short> t1: match->second) {
						int   repeat1, repeatPos1;
						short strand1;
						tie(repeat1, repeatPos1, strand1) = t1;
						if (strand1 == Globals::REVERSE) {
							repeatPos1 = _inputRepeats[repeat1].getSize() - repeatPos1 - Globals::KMER;
						}
						//cout << "\tword: " << lines[0].substr(readPos1, Globals::KMER) << ", pos: " << readPos1 << ", repeat #: " << repeat1 << ", repeat pos: " << repeatPos1 << ", repeat strand: " << strand1 << endl;
						for (tuple < int, int, int, short > t2: positions2) {
							int   readPos2, repeat2, repeatPos2;
							short strand2;
							tie(readPos2, repeat2, repeatPos2, strand2) = t2;
							if (strand2 == Globals::REVERSE) {
								repeatPos2 = _inputRepeats[repeat2].getSize() - repeatPos2 - Globals::KMER;
							}
							//cout << "\t\tmatching repeats on read2: repeat #: " << repeat2 << ", repeat pos: " << repeatPos2 << ", repeat strand: " << strand2 << ", read pos: " << readPos2 << endl;
							updateCount(repeat1, repeatPos1, readPos1, strand1, lines[0].size(), repeat2, repeatPos2, readPos2, strand2, lines[1].size(), m);
						}
					}
				}
			}
		}
	}
}

void Scaffolder::updateCount (int repeat1, int repeatPos1, int readPos1, short strand1, int lineSize1, int repeat2, int repeatPos2, int readPos2, short strand2, int lineSize2, mutex &m) {
	if (repeat1 == repeat2) {
		return;
	}
	short direction, position;
	int   distance, r1, r2;
	if (strand1 == strand2) {
		//cout << "\t\t\tdirect" << endl;
		direction = Globals::DIRECT;
	}
	else {
		//cout << "\t\t\treverse" << endl;
		direction = Globals::REVERSE;
	}
	if (repeat1 > repeat2) {
		//cout << "\t\t\t1 > 2" << endl;
		r1 = repeat1;
		r2 = repeat2;
		distance = _insertSize - readPos2 + (repeatPos2 - _inputRepeats[repeat2].getSize()) + (readPos1 - (lineSize1 - Globals::KMER)) - (repeatPos1 + Globals::KMER);
		position = Globals::AFTER;
		if (strand1 == Globals::REVERSE) {
			//cout << "\t\t\tinvert pos" << endl;
			position = 1 - position;
		}
	}
	else {
		//cout << "\t\t\t1 < 2" << endl;
		r1 = repeat2;
		r2 = repeat1;
		distance = _insertSize - readPos1 + (repeatPos1 - _inputRepeats[repeat1].getSize()) + (readPos2 - (lineSize2 - Globals::KMER)) - (repeatPos2 + Globals::KMER);
		position = Globals::BEFORE;
		if (strand2 == Globals::REVERSE) {
			//cout << "\t\t\tinvert pos" << endl;
			position = 1 - position;
		}
	}
	//cout << "\t\t\tvalues: repeat #" << r1 << " vs repeat #" << r2 << ", position: " << position << ", direction: " << direction << ", distance: " << distance << endl;
	if (_distances[r1][r2][position][direction].size() < Globals::MAX_SCAFFOLD_COUNTS) {
		lock_guard<mutex> lock(m);
		_distances[r1][r2][position][direction].push_back(distance);
	}
}

/*
void Scaffolder::removeWeakLinks () {
	//for (unsigned int i = 0; i < _distances.size(); i++) {
	for (const auto &distanceIIt: _distances) {
		int i = distanceIIt->first;
		if (! _inputRepeats.isRemoved(i)) {
			//for (unsigned int j = 0; j < _distances[i].size(); j++) {
			for (const auto &distanceIJIt: distanceIIt->second) {
				int j = distanceIJIt->first;
				if (! _inputRepeats.isRemoved(j)) {
					auto distanceIJ = distanceIJIt->second;
					for (int k = 0; k < Globals::POSITIONS; k++) {
						for (int l = 0; l < Globals::DIRECTIONS; l++) {
							if (distanceIJ[k][l].size() < Globals::MIN_SCAFFOLD_KMERS) {
								//cout << "Removing weak link " << _distances[i][j][k][l].size() << " < " << Globals::MIN_SCAFFOLD_KMERS << "." << endl;
								distanceIJ[k][l].clear();
							}
						}
					}
				}
			}
		}
	}
}
*/

int Scaffolder::computeMode (const vector <int> &distances) const {
	if (distances.empty()) return numeric_limits<int>::min();
	map <int, unsigned int> values;
	for (int i: distances) {
		values[i]++;
	}
	unsigned int maxValue = 0, mode = 0;
	for (auto &i: values) {
		//cout << i.first << "x" << i.second << " ";
		if (i.second > maxValue) {
			maxValue = i.second;
			mode     = i.first;
		}
	}
	if (maxValue >= Globals::MIN_SCAFFOLD_KMERS) return mode;
	return numeric_limits<int>::min();
}

/*
void Scaffolder::computeMode () {
	for (unsigned int i = 0; i < _distances.size(); i++) {
		if (! _inputRepeats.isRemoved(i)) {
			for (unsigned int j = 0; j < _distances[i].size(); j++) {
				if (! _inputRepeats.isRemoved(j)) {
					for (int k = 0; k < Globals::POSITIONS; k++) {
						for (int l = 0; l < Globals::DIRECTIONS; l++) {
							const vector<int> &v = _distances[i][j][k][l];
							if (! v.empty()) {
								map <int, unsigned int> values;
								//cout << "Values for i: " << i << ", j: " << j << ", position: " << position << ", direction: " << direction << endl;
								for (int i: v) {
									values[i]++;
								}
								unsigned int maxValue = 0, mode = 0;
								for (auto &i: values) {
									//cout << i.first << "x" << i.second << " ";
									if (i.second > maxValue) {
										maxValue = i.second;
										mode     = i.first;
									}
								}
								//cout << endl;
								if (maxValue >= Globals::MIN_SCAFFOLD_KMERS) {
									_modes[i][j][k][l] = mode;
									//cout << "Distance is fine: " << static_cast<int>(static_cast<float>(accumulate(v.begin(), v.end(), 0)) / static_cast<float>(v.size())) << endl;
									//return static_cast<int>(static_cast<float>(accumulate(v.begin(), v.end(), 0)) / static_cast<float>(v.size()));
								}
								//cout << "Modes are " << val1 << " with " << mode1 << " vs " << val2 << " with " << mode2 << endl;
							}
						}
					}
				}
			}
		}
	}
}
*/

/*
void Scaffolder::checkDistances () {
	for (unsigned int i = 0; i < _distances.size(); i++) {
		if (! _inputRepeats.isRemoved(i)) {
			for (unsigned int j = 0; j < _distances[i].size(); j++) {
				if (! _inputRepeats.isRemoved(j)) {
					for (int k = 0; k < Globals::POSITIONS; k++) {
						for (int l = 0; l < Globals::DIRECTIONS; l++) {
							int distance = _modes[i][j][k][l];
							if (distance != numeric_limits<int>::min()) {
								//int distance = static_cast<int>(static_cast<float>(accumulate(v.begin(), v.end(), 0)) / static_cast<float>(v.size()));
								if ((distance < 0) && ((_inputRepeats[i].getSize() < static_cast<unsigned int>(-distance)) || (_inputRepeats[j].getSize() < static_cast<unsigned int>(-distance)))) {
									_modes[i][j][k][l] = numeric_limits<int>::min();
									_distances[i][j][k][l].clear();
									//cout << "\tDistance is too small: repeat1 is " << _inputRepeats[i].getSize() << ", " << "repeat2 is " << _inputRepeats[j].getSize() << ", distance is " << distance << endl;
								}
							}
						}
					}
				}
			}
		}
	}
}
*/

void Scaffolder::resetMinEvidences () {
	unsigned int minEvidences = 1;
	unsigned int maxEvidences = 1;
	unsigned int middleEvidences = 1;
	//for (unsigned int i = 0; i < _distances.size(); i++) {
	for (const auto &distanceIIt: _distances) {
		int i = distanceIIt.first;
		if (! _inputRepeats.isRemoved(i)) {
			//for (unsigned int j = 0; j < _distances[i].size(); j++) {
			for (const auto &distanceIJIt: distanceIIt.second) {
				int j = distanceIJIt.first;
				if (! _inputRepeats.isRemoved(j)) {
					auto &distanceIJ = distanceIJIt.second;
					for (int k = 0; k < Globals::POSITIONS; k++) {
						for (int l = 0; l < Globals::DIRECTIONS; l++) {
							maxEvidences = max<unsigned int>(maxEvidences, distanceIJ[k][l].size());
						}
					}
				}
			}
		}
	}
	if (! aboveEvidenceThreshold(minEvidences)) {
		return;
	}
	middleEvidences = (minEvidences + maxEvidences) / 2;
	while (minEvidences < maxEvidences) {
		if (aboveEvidenceThreshold(middleEvidences)) {
			minEvidences = middleEvidences+1;
		}
		else {
			maxEvidences = middleEvidences-1;
		}
		middleEvidences = (minEvidences + maxEvidences) / 2;
	}
	_maxEvidencesPerNode = middleEvidences;
}

void Scaffolder::buildGraphs() {
	set    < unsigned int > pastNodes;
	vector < unsigned int > currentNodes;
	for (unsigned int start = 0; start < _inputRepeats.getNbRepeats(); start++) {
		if ((! _inputRepeats.isRemoved(start)) && (pastNodes.find(start) == pastNodes.end())) {
			currentNodes.clear();
			SequenceGraph graph;
			graph.addNode(0, _inputRepeats[start].getCount(), _inputRepeats[start].getRepeat());
			currentNodes.push_back(start);
			//cout << "\t\tadding first: " << start << endl;
			for (unsigned int currentIndex = 0; currentIndex < currentNodes.size(); currentIndex++) {
				unsigned int i = currentNodes[currentIndex];
				for (unsigned int j = 0; j < _inputRepeats.getNbRepeats(); j++) {
					if ((j != i) && (! _inputRepeats.isRemoved(j))) {
						for (int k = 0; k < Globals::POSITIONS; k++) {
							for (int l = 0; l < Globals::DIRECTIONS; l++) {
								if ((checkCell(i, j, k, l)) && ((getCell(i, j, k, l).size() >= max<unsigned int>(Globals::MIN_SCAFFOLD_KMERS, _maxEvidencesPerNode)))) {
									int mode = computeMode(getCell(i, j, k, l));
									if ((mode != numeric_limits<int>::min()) && ((mode >= 0) || ((_inputRepeats[i].getSize() >= static_cast<unsigned int>(-mode)) && (_inputRepeats[j].getSize() >= static_cast<unsigned int>(-mode))))) {
										//cout << "\t\tadding: " << i << ", " << j << ", " << k << ", " << l << endl;
										//vector<unsigned int>::iterator end = (currentIndex+1 == currentNodes.size())? currentNodes.end(): currentNodes.begin()+currentIndex+1;
										//vector<unsigned int>::iterator it = find(currentNodes.begin(), end, j);
										vector<unsigned int>::iterator it = find(currentNodes.begin(), currentNodes.end(), j);
										unsigned int nextIndex;
										//cout << "data:"; for (unsigned int c: currentNodes) cout << " " << c; cout << endl;
										if (it == currentNodes.end()) {
											//cout << "\t\t\tfor good" << endl;
											nextIndex = currentNodes.size();
											graph.addNode(nextIndex, _inputRepeats[j].getCount(), _inputRepeats[j].getRepeat());
											currentNodes.push_back(j);
										}
										else {
											nextIndex = it - currentNodes.begin();
										}
										graph.addLink(currentIndex, k, l, nextIndex);
									}
								}
							}
						}
					}
				}
			}
			//cout << "Got this graph:\n" << graph << endl;
			findPaths(graph, currentNodes);
			pastNodes.insert(currentNodes.begin(), currentNodes.end());
		}
	}
}

void Scaffolder::findPaths(SequenceGraph &graph, const vector <unsigned int> &translation) {
	//cout << "Before trimming" << endl;
	//cout << graph << endl;
	GraphTrimmer trimmer(graph);
	//trimmer.removeHubs();
	//cout << "After hub removal" << endl;
	trimmer.trim(false);
	//cout << "After trimming" << endl;
	//cout << graph << endl;
	//cout << "Graph set" << endl;
	graph.findAllPathes();
	//cout << "Paths found" << endl;
	graph.selectPathes(0);
	//cout << "Paths selected" << endl;
	//cout << "Graph is:\n" << graph << endl;
	const vector <SequencePath> &pathes = graph.getPathes();
	for (const SequencePath &path: pathes) {
		_outputRepeats.addRepeat(mergePath(graph, path, translation));
	}
}

CountedRepeat Scaffolder::mergePath(SequenceGraph &graph, const SequencePath &path, const vector <unsigned int> &translation) {
	unsigned int        currentId     = path.getNode(0);
	const SequenceNode &currentNode   = graph.getNode(currentId);
	string              currentString = currentNode.getSequence().getFirstWord();
	int                 currentSize   = currentString.size();
	KmerNb              currentCount  = currentNode.getCount();
	short               position      = -1;
	short               direction     = -1;
	//cout << "Current path is " << path << endl;
	//cout << "\tFirst string: " << currentString << endl;
	for (unsigned int i = 1; i < path.getSize(); i++) {
		unsigned int        nextId        = path.getNode(i);
		short               nextPosition  = path.getNodePosition(i);
		short               nextDirection = path.getNodeDirection(i);
		const SequenceNode &nextNode      = graph.getNode(nextId);
		int                 nextSize      = nextNode.getSequence().getSize();
		KmerNb              nextCount     = nextNode.getCount();
		if (position == -1) {
			position  = nextPosition;
			direction = nextDirection;
		}
		else {
			direction = (nextDirection == Globals::DIRECT)? direction: 1-direction;
		}
		string       nextString = nextNode.getSequence().getWord(direction), firstString, secondString, fill;
		int          thisFirstId, thisSecondId;
		short        thisPosition, thisDirection;
		unsigned int translatedCurrentId = translation[currentId], translatedNextId = translation[nextId];
		if (translatedCurrentId > translatedNextId) {
			thisFirstId   = translatedCurrentId;
			thisSecondId  = translatedNextId;
			thisPosition  = nextPosition;
			thisDirection = nextDirection;
		}
		else {
			thisFirstId   = translatedNextId;
			thisSecondId  = translatedCurrentId;
			thisPosition  = (nextDirection == Globals::DIRECT)? 1-nextPosition: nextPosition;
			thisDirection = nextDirection;
		}
		firstString  = (position == Globals::AFTER)? currentString: nextString;
		secondString = (position == Globals::AFTER)? nextString: currentString;
		/*
		firstString  = ((currentId > nextId) == (position == Globals::AFTER) == (nextDirection == Globals::DIRECT))? currentString: nextString;
		secondString = ((currentId > nextId) == (position == Globals::AFTER) == (nextDirection == Globals::DIRECT))? nextString: currentString;
		*/
		int distance = computeMode(getCell(thisFirstId, thisSecondId, thisPosition, thisDirection));
		//int distance = _modes[thisFirstId][thisSecondId][thisPosition][thisDirection];
		//cout << "distance 1: " << distance << endl;
		if (distance < 0) {
			//cout << "distance: " << distance << endl;
			int accurateDistance = stitch(firstString, secondString, _distances[thisFirstId][thisSecondId][thisPosition][thisDirection]);
			//cout << "Stitched. Check that " << firstString.substr(firstString.size()-accurateDistance, accurateDistance) << " is somehow similar to " << secondString.substr(0, accurateDistance) << endl;
			if (accurateDistance > 0) {
				secondString = secondString.substr(accurateDistance);
			}
			distance = 0;
		}
		//cout << "distance 2: " << distance << endl;
		if (distance > 0) {
			fill = string(distance, 'N');
		}
		currentString = firstString + fill + secondString;
		currentCount  = (currentCount * currentSize + nextCount * nextSize) / (currentSize + nextSize);
		currentSize  += nextSize;
		currentId     = nextId;
		//cout << "\tAdding strings: " << firstString << ", " << secondString << ", currentId: " << currentId << ", nextId: " << nextId << ", thisPosition: " << thisPosition << ", thisDirection " << thisDirection << ", nextPosition: " << nextPosition << ", nextDirection: " << nextDirection << endl;
	}
	//cout << "current string: " << currentString << endl;
	return CountedRepeat(currentString, currentCount, path.isCycle());
}

/*
tuple <int, int, short, short> Scaffolder::findBestScaffold () const {
	int smallest = numeric_limits<int>::max();
	int bestI = -1, bestJ, bestPosition, bestDirection;
	for (unsigned int i = 0; i < _distances.size(); i++) {
		if (! _inputRepeats.isRemoved(i)) {
			for (unsigned int j = 0; j < _distances[i].size(); j++) {
				if (! _inputRepeats.isRemoved(j)) {
					for (int k = 0; k < Globals::POSITIONS; k++) {
						for (int l = 0; l < Globals::DIRECTIONS; l++) {
							if (_distances[i][j][k][l].size() > 0) {
								//cout << "i: " << i << ", j: " << j << ", k: " << k << ", l: " << l << ", values: " << _distances[i][j][k][l] << ", " << _distances[i][j][k][l].size() << endl;
								int avgDistance = static_cast<int>(static_cast<float>(accumulate(_distances[i][j][k][l].begin(), _distances[i][j][k][l].end(), 0)) / static_cast<float>(_distances[i][j][k][l].size()));
								if (avgDistance < smallest) {
									smallest      = avgDistance;
									bestI         = i;
									bestJ         = j;
									bestPosition  = k;
									bestDirection = l;
								}
							}
						}
					}
				}
			}
		}
	}
	return make_tuple(bestI, bestJ, bestPosition, bestDirection);
}
*/

/*
bool Scaffolder::mergeSequences (const int i, const int j, const short position, const short direction, int distance) {
	string repeat1;
	string repeat2;
	string fill = "";
	//cout << "merge sequence with i: " << i << ", j: " << j << ", position: " << position << ", direction: " << direction << ", distance: " << distance << endl;
	if (position == Globals::AFTER) {
		//cout << "pos after" << endl;
		repeat1 = _inputRepeats.getRepeat(i).getRepeat().getFirstWord();
		repeat2 = _inputRepeats.getRepeat(j).getRepeat().getWord(direction);
	}
	else {
		//cout << "pos before" << endl;
		repeat1 = _inputRepeats.getRepeat(j).getRepeat().getWord(direction);
		repeat2 = _inputRepeats.getRepeat(i).getRepeat().getFirstWord();
	}
	//cout << "Sequences are '" << repeat1 << "' and '" << repeat2 << "'" << endl;
	if (distance < 0) {
		if (*max_element(_distances[i][j][position][direction].begin(), _distances[i][j][position][direction].end()) < 0) {
			int accurateDistance = stitch(repeat1, repeat2, _distances[i][j][position][direction]);
			//cout << "Stitched. Check that " << repeat1.substr(repeat1.size()-accurateDistance, accurateDistance) << " is somehow similar to " << repeat2.substr(0, accurateDistance) << endl;
			if (accurateDistance < 0) {
				distance = 0;
			}
			else {
				repeat2 = repeat2.substr(accurateDistance);
			}
		}
		else {
			//cout << "Cannot stitch. Check that " << repeat1.substr(repeat1.size()+distance, -distance) << " is somehow similar to " << repeat2.substr(0, -distance) << endl;
			repeat2 = repeat2.substr(-distance);
		}
	}
	if (distance >= 0) {
		fill = string(distance, 'N');
	}
	repeat1 = repeat1 + fill + repeat2;
	//cout << "\t--> " << repeat1 << endl;
	Sequence newSequence(repeat1);
	KmerNb   newCount = (_inputRepeats.getRepeat(i).getCount() * _inputRepeats.getRepeat(i).getSize() + _inputRepeats.getRepeat(j).getCount() * _inputRepeats.getRepeat(j).getSize()) / (_inputRepeats.getRepeat(i).getSize() + _inputRepeats.getRepeat(j).getSize());
	_inputRepeats.removeRepeat(j);
	_inputRepeats.changeRepeat(i, newSequence, newCount);
	if (newSequence.getSize() >= Globals::MAX_TE_SIZE) {
        _inputRepeats.removeRepeat(i);
    }
	return (newSequence.getFirstWord() != repeat1);
}
*/

bool Scaffolder::checkCell(int i, int j, short position, short direction) const {
	if (i < j) {
		swap(i, j);
		if (direction == Globals::DIRECT) position = 1-position;
	}
	const auto &distanceIIt = _distances.find(i);
	if (distanceIIt == _distances.end()) return false;
	const auto &distanceI = distanceIIt->second;
	const auto &distanceIJIt = distanceI.find(j);
	if (distanceIJIt == distanceI.end()) return false;
	const auto &distanceIJ = distanceIJIt->second;
	return (! distanceIJ[position][direction].empty());
}

const vector <int> &Scaffolder::getCell(int i, int j, short position, short direction) {
	if (i > j) {
		return _distances[i][j][position][direction];
	}
	short newPosition = (direction == Globals::DIRECT)? 1-position: position;
	return _distances[j][i][newPosition][direction];
}

void Scaffolder::setCell(const int i, const int j, const short position, const short direction, const vector <int>&values) {
	if (i > j) {
		//_distances[i][j][position][direction] = values;
		_distances[i][j][position][direction].assign(values.begin(), values.end());
	}
	else if (j > i) {
		short newPosition = (direction == Globals::DIRECT)? 1-position: position;
		//_distances[j][i][newPosition][direction]  = values;
		_distances[j][i][newPosition][direction].assign(values.begin(), values.end());
	}
}

void Scaffolder::freeCell(const int i, const int j, const short position, const short direction) {
	if (i > j) {
		_distances[i][j][position][direction].clear();
	}
	else if (j > i) {
		short newPosition = (direction == Globals::DIRECT)? 1-position: position;
		_distances[j][i][newPosition][direction].clear();
	}
}

int Scaffolder::stitch(const string &seq1, const string &seq2, const vector <int> &v) const {
	//cout << "mode: " << ", min: " << (*min_element(v.begin(), v.end())) << ", max: " << (*max_element(v.begin(), v.end())) << endl;
	int minSize = min<int>(0, -(*max_element(v.begin(), v.end())));
	int maxSize = max<int>(0, -(*min_element(v.begin(), v.end())));
	string s1 = seq1.substr(seq1.size()-min<int>(maxSize, seq1.size()), min<int>(maxSize, seq1.size()));
	string s2 = seq2.substr(0, min<int>(maxSize, seq2.size()));
	int maxStart = s1.size() - 1 - (maxSize - minSize);
	int maxEnd   = maxSize - minSize;
	pair <int, int> p = localAlignment(s1, s2, maxStart, maxEnd);
	return p.second;
}

pair <int, int> Scaffolder::localAlignment(const string &s1, const string &s2, const unsigned int maxStart, const unsigned int maxEnd) const {
	unsigned int maxPenalty = s1.size() + s2.size() + 1;
	vector < vector <unsigned int> > table(s1.size()+1);
	for (unsigned int i = 0; i <= s1.size(); i++) {
		table[i] = vector <unsigned int>(s2.size()+1);
		table[i][0] = (i <= maxStart)? 0: maxPenalty;
	}
	for (unsigned int j = 1; j <= s2.size(); j++) {
		table[0][j] = j;
	}
	for (unsigned int i = 1; i <= s1.size(); i++) {
		for (unsigned int j = 1; j <= s2.size(); j++) {
			table[i][j] = min<unsigned int>(min<unsigned int>(table[i-1][j-1] + ((s1[i-1]==s2[j-1])? 0: 1), table[i][j-1]+1), table[i-1][j]+1);
		}
	}
	int unsigned traceJ = -1, minScore = maxPenalty;
	for (unsigned int j = maxEnd+1; j <= s2.size(); j++) {
		if (table[s1.size()][j] < minScore) {
			minScore = table[s1.size()][j];
			traceJ   = j;
		}
	}
	if (minScore == maxPenalty) {
		//cout << "\tWe have a problem in finding the best alignment between '" << s1 << "' and '" << s2 << "'." << endl;
		return make_pair(-1, -1);
	}
	int i = s1.size(), j = traceJ;
	while (j != 0) {
		if (table[i][j] == table[i-1][j]+1) {
			i--;
		}
		else if (table[i][j] == table[i][j-1]+1) {
			j--;
		}
		else {
			i--;
			j--;
		}
		if (i == 0) {
			j = 0;
		}
	}
	return make_pair(i, traceJ-1);
}

bool Scaffolder::aboveEvidenceThreshold (unsigned int threshold) const {
	unsigned int nbNeighbors;
	for (const auto &distanceIIt: _distances) {
		int i = distanceIIt.first;
		if (! _inputRepeats.isRemoved(i)) {
			nbNeighbors = 0;
			for (const auto &distanceIJIt: distanceIIt.second) {
				int j = distanceIJIt.first;
				if (! _inputRepeats.isRemoved(j)) {
					auto &distanceIJ = distanceIJIt.second;
					for (int k = 0; k < Globals::POSITIONS; k++) {
						for (int l = 0; l < Globals::DIRECTIONS; l++) {
							if (distanceIJ[k][l].size() >= threshold) {
								nbNeighbors++;
							}
						}
					}
				}
			}
			if (nbNeighbors > Globals::SCAFFOLD_MAX_EV) {
				return true;
			}
		}
	}
	return false;
}

Repeats &Scaffolder::getRepeats () {
	return _outputRepeats;
}

ostream& operator<<(ostream& output, const Scaffolder& s) {
	output << "Input: " << s._inputRepeats << endl;
	output << "Output: " << s._inputRepeats << endl;
	for (const auto &distanceIIt: s._distances) {
		int i = distanceIIt.first;
		for (const auto &distanceIJIt: distanceIIt.second) {
			int j = distanceIJIt.first;
			const auto &distanceIJ = distanceIJIt.second;
			output << i << " <-> " << j << ":\n";
			for (int k = 0; k < Globals::POSITIONS; k++) {
				for (int l = 0; l < Globals::DIRECTIONS; l++) {
					output << "\t";
					map <int, int>m;
					for (int v: distanceIJ[k][l]) {
						if (m.find(v) == m.end()) {
							m[v] = 1;
						}
						else {
							m[v]++;
						}
					}
					for (auto &v: m) {
						output << v.first << "x" << v.second << " ";
					}
				}
				output << "\n";
			}
		}
	}
	return output;
}
