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

#include <cstdlib>
#include <set>
#include <algorithm>
#include "globals.hpp"
#include "graphRepeatFinder.hpp"
#include "graphTrimmer.hpp"

GraphRepeatFinder::GraphRepeatFinder(SimpleKmerCount &km, const KmerNb threshold): _kmerCount(km), _threshold(threshold) {}

void GraphRepeatFinder::findRepeats () {
	//cout << "Finding repeats..." << endl;
	unsigned int nbSmall         = 0;
	unsigned int initialHashSize = _kmerCount.getSize();
	while (! _kmerCount.empty()) {
		//cout << "Finding most repeated k-mer..." << endl;
		Kmer firstKmer(getFirstKmer());
		//cout << "  done: " << firstKmer << endl;
		//cout << "Building graph..." << endl;
		SequenceGraph graph;
		fillFirstGraph(graph, firstKmer);
		//cout << "First graph: \n" << graph << endl;
		if (graph.isSmall()) {
			if (! Globals::CHECK.empty() && graph.check()) {
				cout << "\t\t\tSequence is in small graph of size " << graph.getSize() << endl;
			}
			//cout << "\t...graph is small" << endl;
			removeKmers();
			nbSmall++;
		}
		else {
			//cout << "\t...graph is not small" << endl;
			//SequenceGraph bestGraph = findBestGraph(firstGraph);
			if (graph.isBig()) {
				if (! Globals::CHECK.empty() && graph.check()) {
					cout << "\t\t\tSequence is in big graph" << endl;
				}
				graph.findGreedyPathes();
			}
			else {
				if (! Globals::CHECK.empty()) {
					graph.check();
				}
				build(graph);
			}
			//decreaseKmers(bestGraph);
			removeKmers();
			addRepeats(graph);
			nbSmall = 0;
		}
		//cout << "  done." << endl;
		//if (i % 10 == 0) {
		//	cout << i << " graphs solved.";
		//	cout << string(30, '\b') << flush;
		//}
		if ((nbSmall > 0) && (nbSmall % 1000 == 0)) {
			cout << "\t\tFound " << nbSmall << " small graphs in a row (" << (_kmerCount.getSize() * 100 / initialHashSize) << "% hash remaining).    ";
			cout << string(80, '\b') << flush;
		}
		if ((Globals::NB_SMALL_GRAPHS != 0) && (nbSmall >= Globals::NB_SMALL_GRAPHS)) {
			cout << "\tFound " << nbSmall << " small graphs in a row (" << (_kmerCount.getSize() * 100 / initialHashSize) << "% hash remaining). Exiting..." << endl;
			break;
		}
	}
	//cout << i << " graphs solved." << endl;
	//gatherRepeats();
}

KmerCode GraphRepeatFinder::getFirstKmer () {
	pair <KmerCode, KmerNb> p;
	do {
		p = _kmerCount.getRandom();
	}
	while ((p.first != Kmer::UNSET) && (p.second < _threshold));
	return p.first;
}

void GraphRepeatFinder::fillFirstGraph (SequenceGraph &graph, const Kmer &firstKmer) {
	vector < int > indices;
	Sequence      firstSequence = firstKmer.getSequence();
	KmerNb        firstCount    = _kmerCount.getCount(firstKmer);
	KmerCode      firstCode     = firstKmer.getFirstCode();
	//cout << "Inserting nodes..." << endl;
	graph.addNode(0, firstCount, firstSequence);
	_kmers.clear();
	indices.push_back(0);
	_kmers.push_back(firstCode);
	while (! indices.empty()) {
		int currentIndex         = indices.back();
		SequenceNode currentNode = graph.getNode(currentIndex);
		KmerCode currentCode     = _kmers[currentIndex];
		Kmer currentKmer         = Kmer(currentCode);
		KmerNb currentCount      = _kmerCount.getCount(currentKmer);
		indices.pop_back();
		for (short position = 0; position < Globals::POSITIONS; position++) {
			for (short nucleotide = 0; nucleotide < Globals::NB_NUCLEOTIDES; nucleotide++) {
				Kmer nextKmer(currentKmer.getCodeNeighbor(nucleotide, position));
				KmerNb   nextCount    = _kmerCount.getCount(nextKmer);
				KmerCode nextCode     = nextKmer.getFirstCode();
				short    direction    = currentKmer.compare(nextKmer, position);
				auto it = find(_kmers.begin(), _kmers.end(), nextCode);
				if (it != _kmers.end()) {
					graph.addLink(currentIndex, position, direction, it - _kmers.begin());
				}
				else if ((nextKmer.isSet()) && (nextCount >= _threshold) && (nextCount >= currentCount / Globals::FREQUENCY_DIFFERENCE) && (nextCount <= currentCount * Globals::FREQUENCY_DIFFERENCE)) {
				//else if ((nextKmer.isSet()) && (nextCount >= _threshold)) {
					//cout << "Adding " << currentKmer << " <-> " << nextKmer <<  " (position: " << position << ", nucleotide: " << nucleotide << ")" << endl;
					int nextIndex = _kmers.size();
					graph.addNode(nextIndex, nextCount, nextKmer.getSequence());
					indices.push_back(nextIndex);
					graph.addLink(currentIndex, position, direction, nextIndex);
					_kmers.push_back(nextCode);
				}
			}
		}
		if (_kmers.size() % 1000 == 0) {
			cout << "\tBuilding graph with " << _kmers.size() << " nodes explored and " << indices.size() << " in stack.    ";
			cout << string(80, '\b') << flush;
		}
	}
	if (! graph.isSmall()) {
		cout << "\tBuilt graph with " << _kmers.size() << " nodes.                                         " << endl; 
	}
}

/*
Graph GraphRepeatFinder::fillFirstGraph (const Kmer &firstKmer) {
	vector < int > indices;
	Graph        graph(_size, _threshold);
	Sequence     firstSequence = firstKmer.getSequence();
	KmerNb       firstCount    = _kmerCount.getCount(firstKmer);
	int          firstIndex    = graph.addNode(firstSequence, firstCount);
	KmerCode     firstCode     = firstKmer.getFirstCode();
	//cout << "Inserting nodes..." << endl;
	_kmers.clear();
	indices.push_back(firstIndex);
	_kmers.push_back(firstCode);
	while (! indices.empty()) {
		int currentIndex         = indices.back();
		Node currentNode         = graph.getNode(currentIndex);
		Sequence currentSequence = currentNode.getRepeat();
		Kmer currentKmer         = Kmer(currentSequence.getFirstWord());
		KmerNb currentCount      = _kmerCount.getCount(currentKmer);
		indices.pop_back();
		for (short direction = 0; direction < Globals::DIRECTIONS; direction++) {
			for (short nucleotide = 0; nucleotide < Globals::NB_NUCLEOTIDES; nucleotide++) {
				Kmer nextKmer(currentKmer.getCodeNeighbor(nucleotide, direction), _size);
				Sequence nextSequence = nextKmer.getSequence();
				KmerNb   nextCount    = _kmerCount.getCount(nextKmer);
				KmerCode nextCode     = nextKmer.getFirstCode();
				int      found        = -1;
				auto it = find(_kmers.begin(), _kmers.end(), nextCode);
				if (it != _kmers.end()) {
					graph.link(currentIndex, it - _kmers.begin());
				}
				//else if ((nextKmer.isSet()) && (nextCount >= _threshold) && (nextCount >= currentCount / 2)) {
				else if ((nextKmer.isSet()) && (nextCount >= _threshold) && (nextCount >= currentCount / Globals::FREQUENCY_DIFFERENCE) && (nextCount <= currentCount * Globals::FREQUENCY_DIFFERENCE)) {
				//else if ((nextKmer.isSet()) && (nextCount >= _threshold)) {
					//cout << "Adding " << currentKmer << " <-> " << nextKmer <<  " (direction: " << direction << ", nucleotide: " << nucleotide << ")" << endl;
					int nextIndex = graph.addNode(nextSequence, nextCount);
					//_kmerCount.remove(nextKmer);
					indices.push_back(nextIndex);
					graph.link(currentIndex, nextIndex);
					_kmers.push_back(nextCode);
				}
			}
		}
		if (_kmers.size() % 1000 == 0) {
			cout << "\tBuilding graph with " << _kmers.size() << " nodes explored and " << indices.size() << " in stack.";
			cout << string(80, '\b') << flush;
		}
		if (_kmers.size() + indices.size() > Globals::MAX_NB_NODES) {
			cout << "Got a TE of size " << _kmers.size() + indices.size() << "+, which is too much. Please increase minimum depth.\nExiting..." << endl;
			exit(0);
		}
	}
	if (! graph.isSmall()) {
		cout << "\tBuilt graph with " << _kmers.size() << " nodes.                                    " << endl;
	}
	return graph;
}
*/

/*
Graph GraphRepeatFinder::findBestGraph(Graph &firstGraph) {
	firstGraph.collapse();
	//return firstGraph;
	//cout << "First graph:\n" << firstGraph;
	if (firstGraph.getNbBubbles() <= Globals::MAX_BUBBLES) {
		//cout << "\tthis graph seems fine with " << firstGraph.getNbBubbles() << " bubbles" << endl;
		return firstGraph;
	}
	//cout << "\tthis graph has " << firstGraph.getNbBubbles() << " bubbles" << endl;
	KmerNb bestCount     = -1;
	KmerNb minCount      = _threshold;
	KmerNb maxCount      = firstGraph.getMaxCount();
	Graph  lastGoodGraph = firstGraph;
	while (minCount <= maxCount) {
		KmerNb currentCount  = (maxCount + minCount) / 2;
		//cout << "\tbuilding graph with threshold " << currentCount << " between " << minCount << " and " << maxCount << endl;
		Graph currentGraph(firstGraph);
		currentGraph.setThreshold(currentCount);
		currentGraph.collapse();
		//currentGraph.build();
		bool small = currentGraph.isSmall();
		//cout << "\t\t" << currentGraph.getNbBubbles() << " bubbles found. Graph is " << (small? "": "not ") << "small." << endl;
		if ((currentGraph.getNbBubbles() > Globals::MAX_BUBBLES) && (! small)) {
			//cout << "\t\t\tgoing up" << endl;
			minCount = currentCount + 1;
		}
		else {
			//cout << "\t\t\tgoing down" << endl;
			maxCount = currentCount - 1;
		}
		if ((! small) && (currentCount < bestCount) && (currentGraph.getNbBubbles() <= Globals::MAX_BUBBLES)) {
			lastGoodGraph = currentGraph;
			bestCount     = currentCount;
		}
	}
	//cout << "\tfound best graph with threshold " << bestCount << endl;
	if (bestCount == -1) {
		return firstGraph;
	}
	return lastGoodGraph;
}
*/

void GraphRepeatFinder::build(SequenceGraph &graph) {
	GraphTrimmer trimmer(graph);
	trimmer.trim(true);
	//cout << "Graph: " << graph << endl;
	if (graph.findAllPathes()) {
		solveProblem(graph);
		graph.selectPathes(_threshold);
	}
}

void GraphRepeatFinder::removeKmers() {
	for (KmerCode &code: _kmers) {
		_kmerCount.remove(code);
	}
}

/*
void GraphRepeatFinder::decreaseKmers(Graph &graph) {
	vector < pair <int, KmerNb> > kmerIds = graph.getUsedKmers();
	for (auto &kmerId: kmerIds) {
		_kmerCount.decreaseNb(_kmers[kmerId.first], kmerId.second);
	}
}
*/

void GraphRepeatFinder::solveProblem(SequenceGraph &graph) {
	vector <SequencePath> &paths = graph.getPathes();
	set <unsigned int> nodeIds;
	EquationSystem equations(paths.size());
	for (unsigned int i = 0; i < paths.size(); i++) {
		for (unsigned int j = 0; j < paths[i].getSize(); j++) {
			nodeIds.insert(graph.getNode(paths[i].getNode(j)).getId());
		}
	}
	for (unsigned int nodeId: nodeIds) {
		equations.addNode(nodeId, graph.getNode(nodeId).getCount());
	}
	for (unsigned int i = 0; i < paths.size(); i++) {
		for (unsigned int j = 0; j < paths[i].getSize(); j++) {
			equations.addNodePath(i, graph.getNode(paths[i].getNode(j)).getId());
		}
	}
	equations.solveSimplex();
	for (unsigned int i = 0; i < paths.size(); i++) {
		paths[i].setCount(equations.getValue(i));
	}
}

/*
void GraphRepeatFinder::gatherRepeats() {
	for (Graph &graph: _graphs) {
		graph.buildRepeats();
		_repeats.addRepeats(graph.getRepeats());
	}
}
*/

void GraphRepeatFinder::addRepeats(SequenceGraph &graph) {
	for (const CountedRepeat &repeat: graph.getRepeats()) {
		_repeats.addRepeat(repeat);
	}
}

Repeats &GraphRepeatFinder::getRepeats() {
	_repeats.cutToThreshold(_threshold);
	_repeats.sort();
	return _repeats;
}
