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

#include <set>
#include <algorithm>
#include <mutex>
#include <thread>
#include "globals.hpp"
#include "repeatMerger.hpp"
#include "graphTrimmer.hpp"

RepeatMerger::RepeatMerger(const Repeats &repeats, const KmerNb minCount): _size(repeats.getNbRepeats()), _minCount(minCount), _inputRepeats(repeats), _maxPenalty(numeric_limits<Penalty>::max()), _comparisons(_size)  {}

void RepeatMerger::addRepeat (const Sequence &repeat, const KmerNb count) {
	_inputRepeats.addRepeat(repeat, count);
}

void RepeatMerger::mergeRepeats () {
	cout << "Merging repeats (" << _size << " elements)..." << endl;
	for (unsigned int i = 0; i < _size; i++) {
		if (_inputRepeats.getRepeat(i).getSize() <= Globals::MIN_MERGE_SIZE) {
			_inputRepeats.removeRepeat(i);
		}
	}
	if (_inputRepeats.getNbRepeats() < _size) {
		_size = _inputRepeats.getNbRepeats();
		cout << "\tShrinking to " << _size << " elements (the rest has size < " << Globals::MIN_MERGE_SIZE << ")." << endl;
	}
	_inputRepeats.cutToThreshold(_minCount);
	if (_inputRepeats.getNbRepeats() < _size) {
		_size = _inputRepeats.getNbRepeats();
		cout << "\tShrinking to " << _size << " elements (the rest has frequency < " << _minCount << ")." << endl;
	}
	if (Globals::MERGE_MAX_NB > 0) {
		_inputRepeats.cutToNLongest(Globals::MERGE_MAX_NB);
		if (_inputRepeats.getNbRepeats() < _size) {
			_size = _inputRepeats.getNbRepeats();
			cout << "\tShrinking to " << _size << " elements (cannot handle too many elements)." << endl;
		}
	}
	buildStructure();
	//cout << "\tcomparing all repeats..." << endl;
	fillStructure();
	//cout << "\t... done." << endl;
	//cout << *this << endl;
	if (Globals::MERGE_MAX_NODES != 0) {
		resetMaxPenalty();
		if (_maxPenalty != numeric_limits<Penalty>::max()) {
			cout << "\tMax penalty adjusted to " << _maxPenalty << endl;
		}
	}
	buildGraphs();
	cout << "\tGraphs built." << endl;
	cleanStructure();
	cout << "\tGraphs cleaned." << endl;
	_inputRepeats.sort();
	cout << "\tdone." << endl;
	//cout << "\t" << cpt << " merges." << endl;
	//cout << *this << endl;
}

void RepeatMerger::buildStructure () {
	_kmerSets = new KmerSet[_size];
	for (unsigned int i = 0; i < _size; i++) {
		_kmerSets[i].setSequence(_inputRepeats[i].getRepeat().getFirstWord());
	}
}

void RepeatMerger::cleanStructure () {
	delete[] _kmerSets;
}

void RepeatMerger::fillStructure () {
	/*
	if (_j != static_cast<unsigned int>(-1)) {
		for (unsigned int i = 0; i < _size; i++) {
			for (int k = 0; k < Globals::DIRECTIONS; k++) {
				for (int l = 0; l < Globals::POSITIONS; l++) {
					if (_j < i) {
						_comparisons[i][_j][k][l].unset();
					}
					else if (_j > i) {
						_comparisons[_j][i][k][l].unset();
					}
				}
			}
		}
	}
	*/
	vector <thread> threads(Globals::NB_THREADS);
	unsigned int i = 1, j = 0, cpt = 0;
	mutex m;
	for (int threadId = 0; threadId < Globals::NB_THREADS; threadId++) {
		threads[threadId] = thread([this, &i, &j, &cpt, &m]() {
			SequenceComparator comparator;
			unsigned int thisI, thisJ, thisCpt, nbRepeats = _size, nComparisons = _size * (_size) / 2;
			while (true) {
				{
					lock_guard<mutex> lock(m);
					thisI = i;
					thisJ = j;
					do {
						j++;
						if (j >= i) {
							i++;
							j = 0;
						}
					}
					while ((i < nbRepeats) && ((_inputRepeats.isRemoved(i)) || (_inputRepeats.isRemoved(j))));
					cpt++;
					thisCpt = cpt;
				}
				if (thisI >= nbRepeats) {
					return;
				}
				compare(comparator, thisI, thisJ);
				if (thisCpt % 10000 == 0) {
					cout << "\t" << thisCpt << "/" << nComparisons << " comparisons." << endl;
				}
			}
		});
	}
	for (int threadId = 0; threadId < Globals::NB_THREADS; threadId++) {
		threads[threadId].join();
	}
	/*
	if (_j == static_cast<unsigned int>(-1)) {
		cout << "\t" << cpt << " comparisons done." << endl;
	}
	*/
}

void RepeatMerger::compare(SequenceComparator &comparator, unsigned int i, unsigned int j) {
	KmerNb countI = _inputRepeats[i].getCount();
	KmerNb countJ = _inputRepeats[j].getCount();
	if (Globals::FREQUENCY_DIFFERENCE * std::min<KmerNb>(countI, countJ) < std::max<KmerNb>(countI, countJ)) {
		return;
	}
	if (! _kmerSets[i].compare(_kmerSets[j])) {
		return;
	}
	string thisWord = _inputRepeats[i].getRepeat().getFirstWord();
	//if (thisWord.size() > Globals::MAX_MERGE_SIZE) {
	//	thisWord = thisWord.substr(thisWord.size() - Globals::MAX_MERGE_SIZE, Globals::MAX_MERGE_SIZE);
	//}
	for (int k = 0; k < Globals::DIRECTIONS; k++) {
		string thatWord = _inputRepeats[j].getRepeat().getWord(k);
		//if (thatWord.size() > Globals::MAX_MERGE_SIZE) {
		//	thatWord = thatWord.substr(0, Globals::MAX_MERGE_SIZE);
		//}
		string thisWord1 = thisWord.substr(max<int>(0, thisWord.size() - Globals::MAX_MERGE_SIZE), Globals::MAX_MERGE_SIZE);
		string thatWord1 = thatWord.substr(0, Globals::MAX_MERGE_SIZE);
		comparator.compare(thisWord1, thatWord1);
		setCell(i, j, k, Globals::AFTER, comparator, thisWord1, thatWord1);
		//cout << thisWord1 << " and " << thatWord1 << " (" << k << ")\n" << comparator << endl;
		string thatWord2 = thatWord.substr(max<int>(0, thatWord.size() - Globals::MAX_MERGE_SIZE), Globals::MAX_MERGE_SIZE);
		string thisWord2 = thisWord.substr(0, Globals::MAX_MERGE_SIZE);
		comparator.compare(thatWord2, thisWord2);
		setCell(i, j, k, Globals::BEFORE, comparator, thatWord2, thisWord2);
		//cout << thatWord2 << " and " << thisWord2 << " (" << k << ")\n" << comparator << endl;
	}
}

void RepeatMerger::resetMaxPenalty () {
	Penalty minPenalty    = numeric_limits<Penalty>::max();
	Penalty maxPenalty    = 0;
	Penalty middlePenalty = 0;
	for (unsigned int i = 1; i < _size; i++) {
		for (unsigned int j = 0; j < i; j++) {
			for (int k = 0; k < Globals::DIRECTIONS; k++) {
				for (int l = 0; l < Globals::POSITIONS; l++) {
					if (isSet(i, j, k, l)) {
						minPenalty = min<Penalty>(maxPenalty, _comparisons[i][j][k][l].getScore());
						maxPenalty = max<Penalty>(maxPenalty, _comparisons[i][j][k][l].getScore());
					}
				}
			}
		}
	}
	if (! underPenaltyThreshold(minPenalty)) {
		return;
	}
	middlePenalty = (minPenalty + maxPenalty) / 2;
	while (minPenalty < maxPenalty) {
		middlePenalty = (minPenalty + maxPenalty) / 2;
		if (underPenaltyThreshold(middlePenalty)) {
			minPenalty = middlePenalty+1;
		}
		else {
			maxPenalty = middlePenalty-1;
		}
		middlePenalty = (minPenalty + maxPenalty) / 2;
	}
	_maxPenalty = middlePenalty;
}

void RepeatMerger::buildGraphs() {
	set    < unsigned int > pastNodes;
	vector < unsigned int > currentNodes;
	for (unsigned int start = 0; start < _size; start++) {
		if ((! _inputRepeats.isRemoved(start)) && (pastNodes.find(start) == pastNodes.end())) {
			//cout << "\t\tStarting with " << start << endl;
			currentNodes.clear();
			SequenceGraph graph;
			graph.addNode(0, _inputRepeats[start].getCount(), _inputRepeats[start].getRepeat());
			currentNodes.push_back(start);
			for (unsigned int currentIndex = 0; currentIndex < currentNodes.size(); currentIndex++) {
				unsigned int i = currentNodes[currentIndex];
				for (unsigned int j = 0; j < _size; j++) {
					if (j != i) {
						for (int k = 0; k < Globals::DIRECTIONS; k++) {
							for (int l = 0; l < Globals::POSITIONS; l++) {
								if ((isSet(i, j, k, l)) && (getComparison(i, j, k, l).getScore() <= _maxPenalty)) {
									//cout << "adding " << i << ", " << j << ", " << k << ", " << l << endl;
									vector<unsigned int>::iterator end = (currentIndex+1 == currentNodes.size())? currentNodes.end(): currentNodes.begin()+currentIndex+1;
									vector<unsigned int>::iterator it = find(currentNodes.begin(), end, j);
									unsigned int nextIndex;
									if (it != end) {
										nextIndex = it - currentNodes.begin();
									}
									else {
										nextIndex = currentNodes.size();
										graph.addNode(nextIndex, _inputRepeats[j].getCount(), _inputRepeats[j].getRepeat());
										currentNodes.push_back(j);
									}
									if (i > j) {
										graph.addLink(currentIndex, l, k, nextIndex);
									}
									else {
										graph.addLink(nextIndex, (k == Globals::DIRECT)? 1-l: l, k, currentIndex);
									}
								}
							}
						}
					}
				}
			}
			findPaths(graph, currentNodes);
			pastNodes.insert(currentNodes.begin(), currentNodes.end());
		}
	}
}

void RepeatMerger::findPaths(SequenceGraph &graph, const vector <unsigned int> &translation) {
	//cout << "Before trimming" << endl;
	//cout << graph << endl;
	GraphTrimmer trimmer(graph);
	//trimmer.removeHubs();
	trimmer.trim(false);
	//cout << "\tAfter trimming" << endl;
	//cout << "\tGraph set" << endl;
	graph.findAllPathes();
	//cout << "\tPaths found" << endl;
	graph.selectPathes(0);
	//cout << "\tPaths selected" << endl;
	//cout << "\tSmall graph is:\n" << graph << endl;
	const vector <SequencePath> &pathes = graph.getPathes();
	for (const SequencePath &path: pathes) {
		_outputRepeats.addRepeat(mergePath(graph, path, translation));
	}
	//cout << "\tAdded repeats." << endl;
}

CountedRepeat RepeatMerger::mergePath(SequenceGraph &graph, const SequencePath &path, const vector <unsigned int> &translation) {
	int                 currentId     = path.getNode(0);
	const SequenceNode &currentNode   = graph.getNode(currentId);
	string              currentString = currentNode.getSequence().getFirstWord();
	int                 currentSize   = currentString.size();
	KmerNb              currentCount  = currentNode.getCount();
	short               position      = -1;
	short               direction     = -1;
	//cout << "\tMerging path " << path << endl;
	//cout << "\t\tstarting with " << currentId << endl;
	for (unsigned int i = 1; i < path.getSize(); i++) {
		int    nextId                = path.getNode(i);
		short  nextPosition          = path.getNodePosition(i);
		short  nextDirection         = path.getNodeDirection(i);
		const SequenceNode &nextNode = graph.getNode(nextId);
		int    nextSize              = nextNode.getSequence().getSize();
		KmerNb nextCount             = nextNode.getCount();
		if (position == -1) {
			position  = nextPosition;
			direction = nextDirection;
		}
		else {
			direction = (nextDirection == Globals::DIRECT)? direction: 1-direction;
		}
		string nextString = nextNode.getSequence().getWord(direction);
		unsigned int thisFirstId, thisSecondId;
		short        thisPosition, thisDirection;
		unsigned int translatedCurrentId = translation[currentId], translatedNextId = translation[nextId];
		//if (currentId > nextId) {
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
		//cout << "\t\tadding " << nextId << ", position " << nextPosition << ", direction " << nextDirection << endl;
		//for (unsigned int element: translation) cout << element << " ";
		//cout << endl;
		//cout << path << endl;
		//cout << "current: " << currentId << ", next: " << nextId << ", first: " << thisFirstId << ", second: " << thisSecondId << ", direction: " << thisDirection << ", position: " << thisPosition << endl;
		//cout << (*this) << endl;
		//cout << "first: " << thisFirstId << ", " << translation[thisFirstId] << ", second: " << thisSecondId << ", " << translation[thisSecondId] << ", " << _comparisons[translation[thisFirstId]] << " " << _comparisons[translation[thisFirstId]][translation[thisSecondId]] << " " << _comparisons[translation[thisFirstId]][translation[thisSecondId]][thisDirection] << " " << _comparisons[translation[thisFirstId]][translation[thisSecondId]][thisDirection][thisPosition] << endl;
		const ComparisonData &data = _comparisons[thisFirstId][thisSecondId][thisDirection][thisPosition];
		//cout << "Got " << data << endl;
		int             startFirst = data.getStartFirst();
		int             endSecond  = data.getEndSecond();
		/*
		int             pos;
		if (translatedCurrentId > translatedNextId) {
			if (currentId > nextId) {
				if (thisPosition == Globals::AFTER) {
					cout << "CASE 1" << endl;
					pos = endSecond;
				}
				else {
					cout << "CASE 2" << endl;
					pos = startFirst;
				}
			}
			else {
				if (thisPosition == Globals::AFTER) {
					cout << "CASE 3" << endl;
					pos = startFirst;
				}
				else {
					cout << "CASE 4" << endl;
					pos = endSecond;
				}
			}
		}
		else {
			if (currentId > nextId) {
				if (thisPosition == Globals::AFTER) {
					cout << "CASE 5" << endl;
					pos = startFirst;
				}
				else {
					cout << "CASE 6" << endl;
					pos = endSecond;
				}
			}
			else {
				if (thisPosition == Globals::AFTER) {
					cout << "CASE 7" << endl;
					pos = endSecond;
				}
				else {
					cout << "CASE 8" << endl;
					pos = startFirst;
				}
			}
		}
		*/
		int             pos        = ((translatedCurrentId > translatedNextId) == (thisPosition == Globals::AFTER))? endSecond: startFirst;
		//int             pos        = ((currentId > nextId) == (thisPosition == Globals::AFTER))? endSecond: startFirst;
		//cout << "With " << nextSize << " and " << pos << endl;
		if (position == Globals::AFTER) {
			//cout << "\t\t\tmerging " << currentString << " and " << nextString << endl;
			if (thisPosition == Globals::AFTER) {
				//cout << "case 1" << endl;
				currentString += nextString.substr(pos, nextSize - pos);
			}
			else {
				//cout << "case 2 with " << nextSize << " and " << pos << endl;
				//cout << "cs: " << currentString << ", ns: " << nextString << ", nz: " << nextSize << ", pos: " << pos << ", sf: " << startFirst << ", es: " << endSecond << endl;
				currentString += nextString.substr(nextSize - pos, pos);
			}
		}
		else {
			//cout << "\t\t\tmerging " << nextString << " and " << currentString << endl;
			if (thisPosition == Globals::BEFORE) {
				//cout << "case 3" << endl;
				currentString = nextString.substr(0, pos) + currentString;
			}
			else {
				//cout << "case 4" << endl;
				currentString = nextString.substr(nextSize - (pos+1), pos+1) + currentString;
			}
		}
		currentCount = (currentCount * currentSize + nextCount * nextSize) / (currentSize + nextSize);
		currentSize += nextSize;
		currentId    = nextId;
	}
	return CountedRepeat(currentString, currentCount, path.isCycle());
}

/*
void RepeatMerger::findBestMerge () {
	Penalty bestScore = Globals::MAX_PENALTY;
	_i = _j = _direction = _position = -1;
	for (unsigned int i = 1; i < _inputRepeats.getNbRepeats(); i++) {
		if (! _inputRepeats[i].empty()) {
			for (unsigned int j = 0; j < i; j++) {
				if ((! _inputRepeats[j].empty())) {
					for (int k = 0; k < Globals::DIRECTIONS; k++) {
						for (int l = 0; l < Globals::POSITIONS; l++) {
							if ((isSet(i, j, k, l)) && (_comparisons[i][j][k][l].getScore() < bestScore)) {
								_i = i;
								_j = j;
								_direction = k;
								_position  = l;
								bestScore = _comparisons[i][j][k][l].getScore();
							}
						}
					}
				}
			}
		}
	}
}

void RepeatMerger::mergeSequences () {
	string thisString = _inputRepeats[_i].getRepeat().getFirstWord();
	string thatString = _inputRepeats[_j].getRepeat().getWord(_direction);
	int endSecond  = _comparisons[_i][_j][_direction][_position].getEndSecond();
	int thisSize = thisString.size();
	int thatSize = thatString.size();
	KmerNb newCount = (_inputRepeats[_i].getCount() * thisSize + _inputRepeats[_j].getCount() * thatSize) / (thisSize + thatSize);
	string firstString  = (_position == Globals::AFTER)? thisString: thatString;
	string secondString = (_position == Globals::AFTER)? thatString: thisString;
	int sizeSecond = (_position == Globals::AFTER)? thatSize: thisSize;
	//cout << "Merging " << firstString << " and " << secondString << " (" << _inputRepeats[_i].getRepeat() << ", " << _inputRepeats[_j].getRepeat() << "), direction: " << _direction << ", position: " << _position << "\t" << _comparisons[_i][_j][_direction][_position] << endl;
	//if (firstString.substr(startFirst, sizeFirst - startFirst) != secondString.substr(0, endSecond)) {
	//	cerr << "Problem while merging " << firstString << " and " << secondString << ", direction: " << _direction << ", position: " << _position << "\nCannot merge " << firstString.substr(startFirst, sizeFirst - startFirst) << " and " << secondString.substr(0, endSecond) << "\n" << _inputRepeats[_i].getRepeat() << "\t" << _inputRepeats[_j].getRepeat() << "\t" << _comparisons[_i][_j][_direction][_position] << endl;
	//}
	secondString = secondString.substr(endSecond, endSecond - sizeSecond);
	firstString += secondString;
	Sequence newSequence(firstString);
	//cout << "  got " << newSequence << endl;
	_inputRepeats.changeRepeat(_i, newSequence, newCount);
	_inputRepeats.removeRepeat(_j);
	if (newSequence.getSize() >= Globals::MAX_TE_SIZE) {
		_inputRepeats.removeRepeat(_i);
	}
	_kmerSets[_i].setSequence(firstString);
}
*/

Repeats &RepeatMerger::getRepeats () {
	return _outputRepeats;
}

bool RepeatMerger::underPenaltyThreshold (unsigned int threshold) const {
	unsigned int nbNeighbors;
	for (unsigned int i = 1; i < _size; i++) {
		if (! _inputRepeats.isRemoved(i)) {
			nbNeighbors = 0;
			for (unsigned int j = 0; j < i; j++) {
				if (! _inputRepeats.isRemoved(j)) {
					for (int k = 0; k < Globals::DIRECTIONS; k++) {
						for (int l = 0; l < Globals::POSITIONS; l++) {
							if ((isSet(i, j, k, l)) && (_comparisons[i][j][k][l].getScore() <= threshold)) {
								nbNeighbors++;
							}
						}
					}
				}
			}
			if (nbNeighbors > Globals::MERGE_MAX_NODES) {
				return true;
			}
		}
	}
	return false;
}

const ComparisonData &RepeatMerger::getComparison(int i, int j, short direction, short position) const {
	if (i > j) {
		return _comparisons[i][j][direction][position];
	}
	else {
		short newPosition = (direction == Globals::DIRECT)? 1-position: position;
		return _comparisons[j][i][direction][newPosition];
	}
}

bool RepeatMerger::isSet(unsigned int i, unsigned int j, short k, short l) const {
	if (i < j) {
		swap(i, j);
		if (k == Globals::DIRECT) {
			l = 1-l;
		}
	}
	if (_comparisons[i].size() <= j) {
		return false;
	}
	return (_comparisons[i][j][k][l].isSet());
}

void RepeatMerger::setCell(unsigned int i, unsigned int j, short k, short l, SequenceComparator &sc, string &first, string &second) {
	if (i < j) {
		swap(i, j);
		if (k == Globals::DIRECT) {
			l = 1-l;
		}
	}
	ComparisonData cd = ComparisonData(sc, first, second);
	if ((cd.getScore() < Globals::MAX_PENALTY) && (cd.getIdentity() >= Globals::MIN_IDENTITY)) {
		if (_comparisons[i].empty()) {
			_comparisons[i].resize(_size - i + 1);
		}
		_comparisons[i][j][k][l] = cd;
	}
}

ostream& operator<<(ostream& output, const RepeatMerger& rm) {
	output << "Input sequences:\n" << rm._inputRepeats;
	output << "Output sequences:\n" << rm._outputRepeats;
	for (unsigned int i = 1; i < rm._size; i++) {
		for (unsigned int j = 0; j < i; j++) {
			for (int k = 0; k < Globals::DIRECTIONS; k++) {
				for (int l = 0; l < Globals::POSITIONS; l++) {
					if (rm.isSet(i, j, k, l)) {
						output << i << ", " << j << ", " << l << ", " << k << ":\t" << rm._comparisons[i][j][k][l] << "\n";
					}
				}
			}
		}
	}
	return output;
}
