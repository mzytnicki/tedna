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
#include "sequenceGraph.hpp"
#include <thread>
#include <stack>
#include <map>

SequenceNode::SequenceNode (unsigned int id, KmerNb count): _id(id), _count(count), _set(true), _marked(false), _direct(true) {}

SequenceNode::SequenceNode (unsigned int id, KmerNb count, const Sequence &sequence): _id(id), _count(count), _sequence(sequence), _set(true), _marked(false), _direct(true) {}

void SequenceNode::addNeighbor(unsigned int id, short position, short direction) {
	if (find(_neighbors[position][direction].begin(), _neighbors[position][direction].end(), id) == _neighbors[position][direction].end()) {
		_neighbors[position][direction].push_back(id);
	}
}

void SequenceNode::setId(unsigned int id) {
	_id = id;
}

unsigned int SequenceNode::getId() const {
	return _id;
}

unsigned int SequenceNode::getSize() const {
	return _sequence.getSize();
}

const Sequence &SequenceNode::getSequence() const {
	return _sequence;
}

void SequenceNode::setSequence(const Sequence sequence) {
	_sequence = sequence;
}

KmerNb SequenceNode::getCount() const {
	return _count;
}

void SequenceNode::unset() {
	_set = false;
}

bool SequenceNode::isSet() const {
	return _set;
}

void SequenceNode::mark() {
	_marked = true;
}

bool SequenceNode::isMarked() const {
	return _marked;
}

bool SequenceNode::isDirect() const {
	return _direct;
}

bool SequenceNode::isLeaf() const {
	return (((_neighbors[Globals::BEFORE][Globals::DIRECT].empty()) && (_neighbors[Globals::BEFORE][Globals::REVERSE].empty())) || ((_neighbors[Globals::AFTER][Globals::DIRECT].empty()) && (_neighbors[Globals::AFTER][Globals::REVERSE].empty())));
}

short SequenceNode::getOnlyDirection() const {
	if ((! _neighbors[Globals::BEFORE][Globals::DIRECT].empty()) || (! _neighbors[Globals::BEFORE][Globals::REVERSE].empty())) {
		return Globals::BEFORE;
	}
	else {
		return Globals::AFTER;
	}
}

int SequenceNode::getNbNeighbors(short position, short direction) const {
	return _neighbors[position][direction].size();
}

const unsigned int SequenceNode::getNeighbor(short position, short direction, int i) const {
	return _neighbors[position][direction][i];
}

bool SequenceNode::merge(short position, short direction, const SequenceNode &n) {
	for (short d = 0; d < Globals::DIRECTIONS; d++) {
		//_neighbors[position][d] = n._neighbors[(direction == Globals::DIRECT)? position: 1-position][(direction == Globals::DIRECT)? d: 1-d];
		const vector <unsigned int> &neighbors = n._neighbors[(direction == Globals::DIRECT)? position: 1-position][(direction == Globals::DIRECT)? d: 1-d];
		_neighbors[position][d].assign(neighbors.begin(), neighbors.end());
	}
	_count = (getCount() * getSize() + n.getCount() * n.getSize()) / (getSize() + n.getSize());
	string thisString = (_direct)? _sequence.getWord(Globals::DIRECT): _sequence.getWord(Globals::REVERSE);
	//string thatString = (_direct == (direction == Globals::DIRECT))? n._sequence.getWord(Globals::DIRECT): n._sequence.getWord(Globals::REVERSE);
	string thatString = ((n._direct) == (direction == Globals::DIRECT))? n._sequence.getWord(Globals::DIRECT): n._sequence.getWord(Globals::REVERSE);
	string firstString, secondString;
	if (position == Globals::AFTER) {
		firstString = thisString;
		secondString = thatString;
	}
	else {
		firstString = thatString;
		secondString = thisString;
	}
	secondString = secondString.substr(Globals::KMER-1);
	string sequence = firstString + secondString;
	_sequence = Sequence(sequence);
	//cout << "\t\tChecking reverse: " << _sequence.getFirstWord() << " vs " << sequence << endl;
	_direct = (_sequence.getFirstWord() == sequence);
	/*
	if (_sequence.getFirstWord() != sequence) {
		_direct = ! _direct;
		return false;
	}
	*/
	return true;
}

void SequenceNode::updateLinks(unsigned int oldId, unsigned int newId, bool direct) {
	for (short p = 0; p < Globals::POSITIONS; p++) {
		for (short d = 0; d < Globals::DIRECTIONS; d++) {
			for (unsigned int n = 0; n < _neighbors[p][d].size(); n++) {
				if (_neighbors[p][d][n] == oldId) {
					if (direct) {
						_neighbors[p][d][n] = newId;
					}
					else {
						_neighbors[p][1-d].push_back(newId);
					}
				}
			}
		}
	}
}

void SequenceNode::removeLinks(short position, short direction) {
	if (direction == Globals::DIRECTIONS) {
		_neighbors[position][Globals::DIRECT].clear();
		_neighbors[position][Globals::REVERSE].clear();
	}
	else {
		_neighbors[position][direction].clear();
	}
}

bool operator==(const SequenceNode& n1, const SequenceNode& n2) {
	return n1.getId() == n2.getId();
}

ostream& operator<<(ostream& output, const SequenceNode& n) {
	if (! n._neighbors[Globals::BEFORE][Globals::DIRECT].empty()) {
		output << "(+) ";
		for (const int neighbor: n._neighbors[Globals::BEFORE][Globals::DIRECT]) {
			output << neighbor << " ";
		}
	}
	if (! n._neighbors[Globals::BEFORE][Globals::REVERSE].empty()) {
		output << "(-) ";
		for (const int neighbor: n._neighbors[Globals::BEFORE][Globals::REVERSE]) {
			output << neighbor << " ";
		}
	}
	output << "<-- (" << n._id << ", " << n._sequence << ", " << n._count << ")";
	if (! n._direct) {
		output << " (-)";
	}
	output << " --> ";
	if (! n._neighbors[Globals::AFTER][Globals::DIRECT].empty()) {
		output << "(+) ";
		for (const int neighbor: n._neighbors[Globals::AFTER][Globals::DIRECT]) {
			output << neighbor << " ";
		}
	}
	if (! n._neighbors[Globals::AFTER][Globals::REVERSE].empty()) {
		output << "(-) ";
		for (const int neighbor: n._neighbors[Globals::AFTER][Globals::REVERSE]) {
			output << neighbor << " ";
		}
	}
	return output;
}


SequencePath::SequencePath (): _cycle(false), _count(-1) {}

SequencePath::SequencePath (const SequencePath &p): _nodes(p._nodes), _nodeIds(p._nodeIds), _cycle(p._cycle), _count(p._count)  {}

SequencePath::SequencePath (unsigned int id): _cycle(false), _count(-1) {
	addNode(id);
}

SequencePath::SequencePath(unsigned int id1, short position, short direction, unsigned int id2): _cycle(false), _count(-1) {
	addNode(id1);
	addNode(position, direction, id2);
}

void SequencePath::addNode(unsigned int id) {
	addNode(-1, -1, id);
}

void SequencePath::addNode(short position, short direction, unsigned int id) {
	_nodes.push_back(make_tuple(position, direction, id));
	_nodeIds.insert(id);
}

void SequencePath::removeLastNode() {
	int id = SequencePath::getLastNode();
	_nodes.pop_back();
	_nodeIds.erase(id);
}

unsigned int SequencePath::getSize() const {
	return _nodes.size();
}

unsigned int SequencePath::getNode(int i) const {
	return get<2>(_nodes[i]);
}

int SequencePath::getNodePosition(int i) const {
	return get<0>(_nodes[i]);
}

int SequencePath::getNodeDirection(int i) const {
	return get<1>(_nodes[i]);
}

unsigned int SequencePath::getLastNode() const {
	return get<2>(_nodes[_nodes.size()-1]);
}

bool SequencePath::empty() const {
	return _nodes.empty();
}

void SequencePath::clear() {
	_nodes.clear();
	_nodeIds.clear();
}

void SequencePath::reverse() {
	int size = _nodes.size();
	vector < tuple <short, short, unsigned int> > nodes(size);
	for (int i = 0; i < size; i++) {
		nodes[i] = make_tuple(-1, -1, get<2>(_nodes[size-1-i]));
	}
	for (int i = 1; i < size; i++) {
		get<0>(nodes[i]) = (get<1>(_nodes[size-i]) == Globals::DIRECT)? 1-get<0>(_nodes[size-i]): get<0>(_nodes[size-i]);
		get<1>(nodes[i]) = get<1>(_nodes[size-i]);
	}
	_nodes = nodes;
}

void SequencePath::setCycle() {
	_cycle = true;
}

bool SequencePath::isCycle() const {
	return _cycle;
}

void SequencePath::setCount(const KmerNb count) {
	_count = count;
}

KmerNb SequencePath::getCount() const {
	return _count;
}

bool SequencePath::contains(int id) const {
	return (_nodeIds.find(id) != _nodeIds.end());
}

bool SequencePath::contains(const SequencePath &p) const {
	for (int node: p._nodeIds) {
		if (! contains(node)) {
			return false;
		}
	}
	return true;
}

bool SequencePath::crosses(const SequencePath &p) const {
	for (int node: p._nodeIds) {
		if (contains(node)) {
			return true;
		}
	}
	return false;
}

bool SequencePath::crossesInterior(const SequencePath &p) const {
	for (unsigned int i = 1; i < p._nodeIds.size()-1; i++) {
		if (contains(p.getNode(i))) {
			return true;
		}
	}
	for (unsigned int i = 1; i < _nodeIds.size()-1; i++) {
		if (p.contains(getNode(i))) {
			return true;
		}
	}
	return false;
}

void SequencePath::trimTo(const unsigned int destination) {
	vector < tuple <short, short, unsigned int> > nodes;
	short direct = true;
	unsigned int pos;
	_nodeIds.clear();
	for (pos = 0; getNode(pos) != destination; pos++) {
		if (getNodeDirection(pos) != Globals::DIRECT) {
			direct = ! direct;
		}
	}
	nodes.push_back(make_tuple(0, 0, destination));
	_nodeIds.insert(destination);
	for (pos++; pos < _nodes.size(); pos++) {
		_nodeIds.insert(getNode(pos));
		nodes.push_back(_nodes[pos]);
	}
	_nodes = nodes;
}

bool SequencePath::canMerge(const SequencePath &p) const {
	if (getLastNode() == p.getNode(0)) {
		//cout << "case A" << endl;
		unsigned int   last1 = getSize()-1;
		short position1  = (getNodeDirection(last1) == Globals::DIRECT)? 1-getNodePosition(last1): getNodePosition(last1);
		//short direction1 = getNodeDirection(last1);
		//return (((position1 == Globals::AFTER) == (direction1 == Globals::DIRECT)) == (p.getNodePosition(1)));
		return (position1 != p.getNodePosition(1));
	}
	else if (p.getLastNode() == getNode(0)) {
		//cout << "case B" << endl;
		unsigned int   last2 = p.getSize()-1;
		short position2  = (p.getNodeDirection(last2) == Globals::DIRECT)? 1-p.getNodePosition(last2): p.getNodePosition(last2);
		//short direction2 = p.getNodeDirection(last2);
		//return (((position2 == Globals::AFTER) == (direction2 == Globals::DIRECT)) == (getNodePosition(1)));
		return (position2 != getNodePosition(1));
	}
	else if (getNode(0) == p.getNode(0)) {
		//cout << "case C" << endl;
		return (getNodePosition(1) != p.getNodePosition(1));
	}
	else if (getLastNode() == p.getLastNode()) {
		//cout << "case D" << endl;
		unsigned int   last1 = getSize()-1;
		unsigned int   last2 = p.getSize()-1;
		short position1 = (getNodeDirection(last1) == Globals::DIRECT)? 1-getNodePosition(last1): getNodePosition(last1);
		//short direction1 = getNodeDirection(last1);
		short position2 = (p.getNodeDirection(last2) == Globals::DIRECT)? 1-p.getNodePosition(last2): p.getNodePosition(last2);
		//short direction2 = p.getNodeDirection(last2);
		return (position1 != position2);
		//return (((position1 == Globals::AFTER) == (direction1 == Globals::DIRECT)) != ((position2 == Globals::AFTER) == (direction2 == Globals::DIRECT)));
	}
	return false;
}

bool SequencePath::formCycle(const SequencePath &p) const {
	if ((getNode(0) != p.getNode(0)) || (getLastNode() != p.getLastNode())) {
		return false;
	}
	if (getNodePosition(1) == p.getNodePosition(1)) {
		return false;
	}
	unsigned int last1 = getSize()-1;
	unsigned int last2 = p.getSize()-1;
	short position1 = (getNodeDirection(last1) == Globals::DIRECT)? 1-getNodePosition(last1): getNodePosition(last1);
	short position2 = (p.getNodeDirection(last2) == Globals::DIRECT)? 1-p.getNodePosition(last2): p.getNodePosition(last2);
	return (position1 != position2);
}

vector < SequencePath > SequencePath::addNode(const SequenceNode &linkNode, const SequenceNode &newNode) const {
	vector < SequencePath > paths;
	unsigned int linkNodeId = linkNode.getId(), newNodeId = newNode.getId();
	if (getSize() == 1) {
		for (short position = 0; position < Globals::POSITIONS; position++) {
			for (short direction = 0; direction < Globals::DIRECTIONS; direction++) {
				for (int i = 0; i < linkNode.getNbNeighbors(position, direction); i++) {
					if (linkNode.getNeighbor(position, direction, i) == newNodeId) {
						SequencePath newPath(linkNodeId);
						newPath.addNode(position, direction, newNodeId);
						paths.push_back(newPath);
					}
				}
			}
		}
	}
	else {
		if (linkNodeId == getNode(0)) {
			short p = getNodePosition(1);
			for (short position = 0; position < Globals::POSITIONS; position++) {
				short direction = Globals::DIRECT;
				if (p != position) {
					direction = Globals::REVERSE;
				}
				for (int i = 0; i < newNode.getNbNeighbors(position, direction); i++) {
					if (newNode.getNeighbor(position, direction, i) == linkNode) {
						SequencePath newPath(newNodeId);
						newPath.addNode(position, direction, linkNodeId);
						for (unsigned int j = 1; j < getSize(); j++) {
							newPath.addNode(getNodePosition(j), getNodeDirection(j), getNode(j));
						}
						paths.push_back(newPath);
					}
				}
			}
		}
		else if (linkNodeId == getLastNode()) {
			unsigned int last = getSize() - 1;
			short p = getNodePosition(last);
			short d = getNodeDirection(last);
			short position = (d == Globals::DIRECT)? p: 1-p;
			for (short direction = 0; direction < Globals::DIRECTIONS; direction++) {
				for (int i = 0; i < linkNode.getNbNeighbors(position, direction); i++) {
					if (linkNode.getNeighbor(position, direction, i) == newNodeId) {
						SequencePath newPath(*this);
						newPath.addNode(position, direction, newNodeId);
						paths.push_back(newPath);
					}
				}
			}
		}
		else {
			cerr << "Problem in 'AddNode'" << endl;
		}
	}
	return paths;
}

SequencePath SequencePath::merge(const SequenceNode &linkNode, const SequencePath &path) const {
	unsigned int linkNodeId = linkNode.getId();
	SequencePath p1 = *this, p2 = path, pOut;
	if ((p1.getNode(0) == linkNodeId) && (p2.getLastNode() == linkNodeId)) {
		swap(p1, p2);
	}
	else {
		if (p1.getNode(0) == linkNodeId) {
			p1.reverse();
		}
		if (p2.getLastNode() == linkNodeId) {
			p2.reverse();
		}
	}
	for (unsigned int nodeId: path._nodeIds) {
		if ((nodeId != linkNodeId) && (_nodeIds.find(nodeId) != _nodeIds.end())) {
			return pOut;
		}
	}
	if ((p1.getLastNode() == linkNodeId) && (p2.getNode(0) == linkNodeId)) {
		unsigned int last1 = p1.getSize()-1;
		short d = p1.getNodeDirection(last1), p = p1.getNodePosition(last1);
		short position1 = (d == Globals::DIRECT)? p: 1-p;
		if (position1 == p2.getNodePosition(1)) {
			for (unsigned int i = 0; i < p1.getSize(); i++) {
				pOut.addNode(p1.getNodePosition(i), p1.getNodeDirection(i), p1.getNode(i));
			}
			for (unsigned int i = 1; i < p2.getSize(); i++) {
				pOut.addNode(p2.getNodePosition(i), p2.getNodeDirection(i), p2.getNode(i));
			}
		}
	}
	else {
		cerr << "Problem while merging paths" << endl;
	}
	return pOut;
}

SequencePath operator+(const SequencePath& p1, const SequencePath& p2) {
	SequencePath pOut, pIn1, pIn2;
	if (p1.getLastNode() == p2.getNode(0)) {
		//cout << " case 1" << endl;
		pIn1 = p1;
		pIn2 = p2;
	}
	else if (p2.getLastNode() == p1.getNode(0)) {
		//cout << " case 2" << endl;
		pIn1 = p2;
		pIn2 = p1;
	}
	else if (p1.getNode(0) == p2.getNode(0)) {
		//cout << " case 3. pIn1: " << pIn1;
		pIn1 = p1;
		pIn2 = p2;
		pIn1.reverse();
		//cout << " then " << pIn1 << endl;
	}
	else if (p1.getLastNode() == p2.getLastNode()) {
		//cout << " case 4" << endl;
		pIn1 = p1;
		pIn2 = p2;
		pIn2.reverse();
	}
	//cout << "+ing " << p1 << " and " << p2 << ", pIn1 is " << pIn1 << " and pIn2 is " << pIn2 << endl;
	for (unsigned int i = 0; i < pIn1.getSize(); i++) {
		pOut.addNode(pIn1.getNodePosition(i), pIn1.getNodeDirection(i), pIn1.getNode(i));
	}
	for (unsigned int i = 1; i < pIn2.getSize(); i++) {
		pOut.addNode(pIn2.getNodePosition(i), pIn2.getNodeDirection(i), pIn2.getNode(i));
	}
	/*
	pOut._nodes.insert(pOut._nodes.end(), p1._nodes.begin(), p1._nodes.end());
	pOut._nodes.insert(pOut._nodes.end(), p2._nodes.begin(), p2._nodes.end());
	for (int i = 0; i < pIn1.getSize(); i++) {
		pOut._nodes.push_back(pIn1._nodes[i]);
	}
	for (int i = 1; i < pIn2.getSize(); i++) {
		pOut._nodes.push_back(pIn2._nodes[i]);
	}
	for (int i: p1._nodeIds) {
		pOut._nodeIds.insert(i);
	}
	for (int i: p2._nodeIds) {
		pOut._nodeIds.insert(i);
	}
	*/
	if (pOut.getNode(0) > pOut.getLastNode()) {
		pOut.reverse();
		//cout << "\t\tReversing" << endl;
	}
	//cout << "\t\tInput paths: " << p1 << " and " << p2 << endl;
	return pOut;
}

bool operator==(const SequencePath& p1, const SequencePath& p2) {
	if (p1.isCycle() != p2.isCycle()) {
		return false;
	}
	if (p1._nodes.size() != p2._nodes.size()) {
		return false;
	}
	return equal(p1._nodeIds.begin(), p1._nodeIds.end(), p2._nodeIds.begin());
}

ostream& operator<<(ostream& output, const SequencePath& p) {
	for (const tuple < short, short, int > &element: p._nodes) {
		int id;
		short position, direction;
		tie(position, direction, id) = element;
		if (position == Globals::AFTER) {
			output << " -->";
		}
		else if (position == Globals::BEFORE) {
			output << " <--";
		}
		if (direction == Globals::DIRECT) {
			output << " (+)";
		}
		if (direction == Globals::REVERSE) {
			output << " (-)";
		}
		output << " " << id;
	}
	if (p.getCount() != static_cast<KmerNb>(-1)) {
		output << " -- " << p.getCount();
	}
	if (p.isCycle()) {
		output << " -- cycle";
	}
	return output;
}


SequenceGraph::SequenceGraph (): _maxPaths(0) {}

SequenceGraph::SequenceGraph (int size): _maxPaths(0), _nodes(size) {}

void SequenceGraph::setMaxPaths(const unsigned int maxPaths) {
	_maxPaths = maxPaths;
}

unsigned int SequenceGraph::getSize() const {
	return _nodes.size();
}

bool SequenceGraph::isSmall () const {
	int size = Globals::KMER-1;
	for (const SequenceNode &node: _nodes) {
		//cout << "Size of node " << node << " is " << node.getSize() << endl;
		if (node.isSet()) {
			//cout << "\tis set" << endl;
			//size += node.getSize();
			size++;
		}
	}
	//cout << "Total size is " << size << ", threshold is " << Globals::MIN_NB_NODES << endl;
	return (size < Globals::MIN_NB_NODES);
}

bool SequenceGraph::isBig () const {
	return (_nodes.size() > Globals::MAX_NB_NODES);
}

void SequenceGraph::addNode(unsigned int id, KmerNb count) {
	if (id < _nodes.size()) {
		_nodes[id] = SequenceNode(id, count);
	}
	else if (id == _nodes.size()) {
		_nodes.push_back(SequenceNode(id, count));
	}
	else {
		cerr << "Error while inserting the node " << id << endl;
	}
}

void SequenceGraph::addNode(unsigned int id, KmerNb count, const Sequence &sequence) {
	if (id < _nodes.size()) {
		//cout << "adding node " << id << "/" << _nodes.size() << " mode 1" << endl;
		_nodes[id] = SequenceNode(id, count, sequence);
	}
	else if (id == _nodes.size()) {
		//cout << "adding node " << id << "/" << _nodes.size() << " mode 2" << endl;
		_nodes.push_back(SequenceNode(id, count, sequence));
	}
	else {
		cerr << "Error while inserting the node " << id << endl;
	}
}

void SequenceGraph::addLink (unsigned int n1, short pos1, short dir1, unsigned int n2) {
	_nodes[n1].addNeighbor(n2, pos1, dir1);
	short dir2 = dir1;
	short pos2 = (dir1 == Globals::DIRECT)? 1-pos1: pos1;
	_nodes[n2].addNeighbor(n1, pos2, dir2);
}

SequenceNode &SequenceGraph::getNode (const unsigned int i) {
	return _nodes[i];
}

bool SequenceGraph::findAllPathes() {
	//findLeaves();
	return findPathes();
	//findCycles();
}

bool SequenceGraph::findPathes() {
	//cout << "Starting path finding with graph\n" << *this << endl;
	map < unsigned int, vector < unsigned int > > pathEnds;
	for (unsigned int nodeId = 0; nodeId < getSize(); nodeId++) {
		SequenceNode &node = _nodes[nodeId];
		if (node.isSet()) {
			//cout << "\tNow node " << node << endl;
			unsigned int start = _pathes.size();
			for (short position = 0; position < Globals::POSITIONS; position++) {
				for (short direction = 0; direction < Globals::DIRECTIONS; direction++) {
					for (int i = 0; i < node.getNbNeighbors(position, direction); i++) {
						unsigned int neighborId = node.getNeighbor(position, direction, i);
						SequenceNode &neighbor = _nodes[neighborId];
						if (neighbor.isSet()) {
							//cout << "\t\tNow neighbor " << neighbor << endl;
							unsigned int peEnd = pathEnds[neighborId].size();
							for (unsigned int pathId = 0; pathId < peEnd; pathId++) {
								SequencePath &path = _pathes[pathEnds[neighborId][pathId]];
								//cout << "\t\t\tNow path " << path << endl;
								for (SequencePath &newPath: path.addNode(neighbor, node)) {
									if (checkPath(newPath)) {
										//cout << "\t\t\t\tNow new path " << newPath << endl;
										if (! newPath.isCycle()) {
											pathEnds[newPath.getNode(0)].push_back(_pathes.size());
											pathEnds[newPath.getLastNode()].push_back(_pathes.size());
										}
										_pathes.push_back(newPath);
										if ((_maxPaths > 0) && (_pathes.size() > _maxPaths)) {
											_pathes.clear();
											findGreedyPathes();
											return false;
										}
									}
								}
							}
						}
					}
				}
			}
			unsigned int end = _pathes.size();
			for (unsigned int pathId = start; pathId < end; pathId++) {
				SequencePath &path = _pathes[pathId];
				if (! path.isCycle()) {
					//cout << "\t\tCycling path " << path << endl;
					unsigned int last = path.getSize() - 1;
					short pStart = path.getNodePosition(1);
					short pEnd   = path.getNodePosition(last);
					short dEnd   = path.getNodeDirection(last);
					short p      = (dEnd == Globals::DIRECT)? pEnd: 1-pEnd;
					short d      = Globals::DIRECT;
					bool found   = false;
					SequenceNode &lastNode   = _nodes[path.getNode(last)];
					unsigned int firstNodeId = path.getNode(0);
					if (p != pStart) {
						d = Globals::REVERSE;
					}
					for (int i = 0; (i < lastNode.getNbNeighbors(p, d)) && (! found); i++) {
						if (lastNode.getNeighbor(p, d, i) == firstNodeId) {
							found = true;
						}
					}
					if (found) {
						SequencePath newPath(path);
						newPath.addNode(p, d, firstNodeId);
						newPath.setCycle();
						if (checkPath(newPath)) {
							//cout << "\t\t\tGot new path " << newPath << endl;
							_pathes.push_back(newPath);
							if ((_maxPaths > 0) && (_pathes.size() > _maxPaths)) {
								_pathes.clear();
								findGreedyPathes();
								return false;
							}
						}
					}
				}
			}
			for (unsigned int pathId1 = start; pathId1 < end; pathId1++) {
				if (! _pathes[pathId1].isCycle()) {
					//cout << "\t\tPath 1 " << _pathes[pathId1] << endl;
					for (unsigned int pathId2 = pathId1+1; pathId2 < end; pathId2++) {
						if (! _pathes[pathId2].isCycle()) {
							//cout << "\t\t\tPath 2 " << _pathes[pathId2] << endl;
							SequencePath newPath = _pathes[pathId1].merge(node, _pathes[pathId2]);
							if ((! newPath.empty()) && (checkPath(newPath))) {
								//cout << "\t\t\t\tNew path " << newPath << endl;
								if (! newPath.isCycle()) {
									pathEnds[newPath.getNode(0)].push_back(_pathes.size());
									pathEnds[newPath.getLastNode()].push_back(_pathes.size());
								}
								_pathes.push_back(newPath);
								if ((_maxPaths > 0) && (_pathes.size() > _maxPaths)) {
									_pathes.clear();
									findGreedyPathes();
									return false;
								}
							}
						}
					}
				}
			}
			pathEnds[nodeId].push_back(_pathes.size());
			_pathes.push_back(SequencePath(nodeId));
			if ((_maxPaths > 0) && (_pathes.size() > _maxPaths)) {
				_pathes.clear();
				findGreedyPathes();
				return false;
			}
		}
	}	
	/*
	cout << "Over with" << endl;
	for (SequencePath &path: _pathes) {
		cout << "\t" << path << endl;
	}
	*/
	//cout << "Found " << _pathes.size() << " pathes." << endl;
	return true;
}

unsigned int SequenceGraph::findMostSeenNode() const {
	KmerNb       highestCount = 0;
	unsigned int highestId    = -1;
	for (unsigned int nodeId = 0; nodeId < getSize(); nodeId++) {
		if ((_nodes[nodeId].isSet()) && (_nodes[nodeId].getCount() > highestCount)) {
			highestCount = _nodes[nodeId].getCount();
			highestId    = nodeId;
		}
	}
	return highestId;
}

pair < unsigned int, bool > SequenceGraph::findMostSeenNeighbor(const unsigned int nodeId, const unsigned short position) const {
	const SequenceNode& node   = _nodes[nodeId];
	KmerNb        highestCount = 0, count;
	unsigned int  highestId    = -1, id;
	short         highestDir   = 0;
	for (short d = 0; d < Globals::DIRECTIONS; d++) {
		for (int i = 0; i < node.getNbNeighbors(position, d); i++) {
			id    = node.getNeighbor(position, d, i);
			count = _nodes[id].getCount();
			if ((_nodes[id].isSet()) && (count > highestCount)) {
				highestCount = count;
				highestId    = id;
				highestDir   = d;
			}
		}
	}
	return make_pair(highestId, highestDir);
}

void SequenceGraph::findGreedyPathes() {
	//cout << "Starting path finding with graph\n" << *this << endl;
	static const unsigned int over = static_cast<unsigned int>(-1);
	unsigned int highestId = findMostSeenNode(), currentId, nextId;
	short        highestDir, p, d;
	while (highestId != over) {
		SequenceNode &highestNode = _nodes[highestId];
		string        sequence    = highestNode.getSequence().getWord(Globals::DIRECT);
		KmerNb        count       = highestNode.getCount();
		unsigned int  cpt         = 1;
		highestNode.unset();
		for (short position = 0; position < Globals::POSITIONS; position++) {
			currentId               = highestId;
			p                       = position;
			d                       = Globals::DIRECT;
			tie(nextId, highestDir) = findMostSeenNeighbor(currentId, p);
			while (nextId != over) {
				currentId = nextId;
				if (highestDir == Globals::REVERSE) {
					p = 1-p;
					d = 1-d;
				}
				SequenceNode &currentNode     = _nodes[currentId];
				string        currentSequence = currentNode.getSequence().getWord(d);
				count                        += currentNode.getCount();
				if (position == Globals::AFTER) {
					sequence += currentSequence.substr(Globals::KMER-1);
				}
				else {
					sequence = currentSequence.substr(0, currentSequence.size() - Globals::KMER + 1) + sequence;
				}
				cpt++;
				currentNode.unset();
				tie(nextId, highestDir) = findMostSeenNeighbor(currentId, p);
			}
		}
		highestId = findMostSeenNode();
		_repeats.push_back(CountedRepeat(sequence, count / cpt));
	}
	/*
	cout << "Over with" << endl;
	for (CountedRepeat &repeat: _repeats) {
		cout << "\t" << repeat << endl;
	}
	*/
}

/*
void SequenceGraph::findPathes() {
	//cout << "Starting path finding" << endl;
	unsigned int nbPaths = _pathes.size()+1;
	while (nbPaths != _pathes.size()) {
		nbPaths = _pathes.size();
		findLeaves();
		vector < unsigned int > leaves;
		vector <thread> threads(Globals::NB_THREADS);
		mutex m1, m2;
		for (unsigned int i = 0; i < _leaves.size(); i++) {
			SequenceNode &node = _nodes[_leaves[i]];
			if (node.isSet()) {
				//cout << "Found leaf " << node << endl;
				leaves.push_back(node.getId());
			}
		}
		for (int threadId = 0; threadId < Globals::NB_THREADS; threadId++) {
			threads[threadId] = thread([this, &leaves, &m1, &m2]() {
				while (true) {
					unsigned int id;
					{
						lock_guard<mutex> lock(m1);
						if (! leaves.empty()) {
							id = leaves.back();
							leaves.pop_back();
						}
						else {
							return;
						}
					}
					findPathes(id, m2);
				}
			});
		}
		for (int threadId = 0; threadId < Globals::NB_THREADS; threadId++) {
			threads[threadId].join();
		}
		removeMarkedNodes();
	}
}
*/

/*
void SequenceGraph::findCycles() {
	mutex m;
	for (const SequenceNode &node: _nodes) {
		if ((node.isSet()) && (! node.isMarked())) {
			//cout << "Starting cycle at " << node << endl;
			findPathes(node.getId(), m);
		}
	}
}
*/

bool SequenceGraph::isLeaf(const SequenceNode &n) const {
	if (n.isLeaf()) {
		return true;
	}
	for (short position = 0; position < Globals::POSITIONS; position++) {
		bool empty = true;
		for (short direction = 0; (direction < Globals::DIRECTIONS) && (empty); direction++) {
			for (int i = 0; (i < n.getNbNeighbors(position, direction)) && (empty); i++) {
				if (_nodes[n.getNeighbor(position, direction, i)].isSet()) {
					empty = false;
				}
			}
		}
		if (empty) {
			return true;
		}
	}
	return false;
}

void SequenceGraph::countNeighbors(const SequenceNode &n, KmerNb *count) const {
	for (short position = 0; position < Globals::POSITIONS; position++) {
		for (short direction = 0; direction < Globals::DIRECTIONS; direction++) {
			for (int i = 0; i < n.getNbNeighbors(position, direction); i++) {
				if (_nodes[n.getNeighbor(position, direction, i)].isSet()) {
					count[position]++;
				}
			}
		}
	}
}

int SequenceGraph::isVee(const SequenceNode &n) const {
	KmerNb count[Globals::POSITIONS] = {0, 0};
	countNeighbors(n, count);
	if ((count[Globals::BEFORE] > 1) && (count[Globals::AFTER] > 1)) {
		return 1+Globals::POSITIONS;
	}
	if (count[Globals::BEFORE] > 1) {
		return 1+Globals::BEFORE;
	}
	if (count[Globals::AFTER] > 1) {
		return 1+Globals::AFTER;
	}
	return 0;
}

int SequenceGraph::isFork(const SequenceNode &n) const {
	if (n.isLeaf()) {
		return false;
	}
	KmerNb count[Globals::POSITIONS] = {0, 0};
	countNeighbors(n, count);
	if ((count[Globals::BEFORE] == 1) && (count[Globals::AFTER] > 1)) {
		return 1 + Globals::AFTER;
	}
	if ((count[Globals::BEFORE] > 1) && (count[Globals::AFTER] == 1)) {
		return 1 + Globals::BEFORE;
	}
	return 0;
}

bool SequenceGraph::isAlone(const SequenceNode &n) const {
	for (short position = 0; position < Globals::POSITIONS; position++) {
		for (short direction = 0; direction < Globals::DIRECTIONS; direction++) {
			for (int i = 0; i < n.getNbNeighbors(position, direction); i++) {
				if (_nodes[n.getNeighbor(position, direction, i)].isSet()) {
					return false;
				}
			}
		}
	}
	return true;
}

short SequenceGraph::getOnlyDirection(const SequenceNode &n) const {
	//cout << "\t\t\t\tGetting directions of " << n << ":";
	for (short direction = 0; direction < Globals::DIRECTIONS; direction++) {
		//cout << " direction: " << direction;
		for (int i = 0; i < n.getNbNeighbors(Globals::BEFORE, direction); i++) {
			//cout << " i: " << i << _nodes[n.getNeighbor(Globals::BEFORE, direction, i)];
			if (_nodes[n.getNeighbor(Globals::BEFORE, direction, i)].isSet()) {
				//cout << endl;
				return Globals::BEFORE;
			}
		}
	}
	//cout << endl;
	return Globals::AFTER;
}

unsigned int SequenceGraph::getCumulatedSize(const SequencePath &path) const {
	unsigned int sum = Globals::KMER-1;
	for (unsigned int i = 0; i < path.getSize(); i++) {
		sum += _nodes[path.getNode(i)].getSize() - (Globals::KMER-1);
	}
	return sum;
}

void SequenceGraph::findLeaves() {
	_leaves.clear();
	for (SequenceNode &node: _nodes) {
		if ((node.isSet()) && (isLeaf(node))) {
			_leaves.push_back(node.getId());
		}
	}
}

void SequenceGraph::removeMarkedNodes() {
	for (SequenceNode &node: _nodes) {
		if (node.isMarked()) {
			node.unset();
		}
	}
}

/*
void SequenceGraph::findPathes() {
	vector < SequencePath > paths;
	vector < vector < unsigned int > > pathEnds(_nodes.size()+1);
	stack < unsigned int > openPaths;
	for (SequenceNode &node: _nodes) {
		if (node.isSet()) {
			int id = node.getId();
			for (short position = 0; position < Globals::POSITIONS; position++) {
				for (short direction = 0; direction < Globals::DIRECTIONS; direction++) {
					for (int neighborId = 0; neighborId < node.getNbNeighbors(position, direction); neighborId++) {
						int neighbor = node.getNeighbor(position, direction, neighborId);
						if ((_nodes[neighbor].isSet()) && (id < neighbor)) {
							SequencePath path(id, position, direction, neighbor);
							pathEnds[id].push_back(paths.size());
							pathEnds[neighbor].push_back(paths.size());
							openPaths.push(paths.size());
							paths.push_back(path);
							//cout << "Adding first " << path << " with size " << path.getSize() << endl;
						}
					}
				}
			}
		}
	}
	while (! openPaths.empty()) {
		unsigned int pathId = openPaths.top();
		SequencePath path = paths[pathId];
		openPaths.pop();
		if (! path.isCycle()) {
			//cout << "\tPicking path1 " << path << " with size " << path.getSize() << endl;
			unsigned int endIds[] = {path.getNode(0), path.getLastNode()};
			for (unsigned int endId: endIds) {
				for (unsigned int otherPathId: pathEnds[endId]) {
					if (otherPathId != pathId) {
						SequencePath otherPath = paths[otherPathId];
						//cout << "\t\tPicking path2 " << otherPath << " with size " << otherPath.getSize() << endl;
						if ((path.canMerge(otherPath)) && (! path.crossesInterior(otherPath))) {
							SequencePath newPath = path + otherPath;
							//cout << "\tSuggesting " << newPath << " with size " << newPath.getSize() << endl;
							//cout << "\t\tStep 1 -- Input pathes are " << path << " and " << otherPath << endl;
							int newStart = newPath.getNode(0);
							int newEnd   = newPath.getLastNode();
							if (newStart == newEnd) {
								if (! path.formCycle(otherPath)) {
									continue;
								}
								newPath.removeLastNode();
								newPath.setCycle();
								newStart = newEnd = _nodes.size();
							}
							//cout << "\t\tStep 2 -- Input pathes are " << path << " and " << otherPath << endl;
							bool seen             = false;
							unsigned int  pathEnd = (pathEnds[newStart].size() < pathEnds[newEnd].size())? newStart: newEnd;
							for (unsigned int storedPathId: pathEnds[pathEnd]) {
								if ((storedPathId != pathId) && (storedPathId != otherPathId)) {
									SequencePath &storedPath = paths[storedPathId];
									if (storedPath.contains(newPath)) {
										//cout << "\tEliminated by " << storedPath << endl;
										seen = true;
										break;
									}
								}
							}
							//cout << "\t\tStep 3 -- Input pathes are " << path << " and " << otherPath << endl;
							if (! seen) {
								pathEnds[newStart].push_back(paths.size());
								if (newStart != newEnd) {
									pathEnds[newEnd].push_back(paths.size());
									openPaths.push(paths.size());
								}
								paths.push_back(newPath);
								//cout << "Adding then " << newPath << " with size " << newPath.getSize() << " and input " << path << " and " << otherPath << endl;
								if (paths.size() % 100 == 0) {
									cout << "\tGot " << paths.size() << " paths (and " << openPaths.size() << " to inspect).";
									cout << string(70, '\b') << flush;
								}
							}
						}
					}
				}
			}
		}
	}
	if (paths.size() >= 1000) {
		cout << "\tGot " << paths.size() << " paths.                            " << endl;
	}
	for (const SequencePath &path: paths) {
		if (((_nodes[path.getNode(0)].isLeaf()) && (_nodes[path.getLastNode()].isLeaf())) || (path.isCycle())) {
			_pathes.push_back(path);
		}
	}
	for (unsigned int i = 0; i < _nodes.size(); i++) {
		if (_nodes[i].isSet()) {
			SequencePath path(i);
			_pathes.push_back(path);
		}
	}
	//cout << "Paths:\n" << endl;
	//for (SequencePath &path: _pathes)
	//	cout << path << endl;
}
*/

/*
void SequenceGraph::findPathes() {
	vector < PathEnds > paths;
	vector < SequencePath > loops;
	Matrix matrix (_nodes.size());
	long nbPaths = 0;
	for (SequenceNode &node: _nodes) {
		int id = node.getId();
		for (short position = 0; position < Globals::POSITIONS; position++) {
			for (short direction = 0; direction < Globals::DIRECTIONS; direction++) {
				for (int neighborId = 0; neighborId < node.getNbNeighbors(position, direction); neighborId++) {
					int neighbor = node.getNeighbor(position, direction, neighborId);
					if (id < neighbor) {
						SequencePath path(id, position, direction, neighbor);
						PathEnds pe(id, neighbor);
						matrix.addPath(id, neighbor, path);
						paths.push_back(pe);
						nbPaths++;
						//cout << "Adding first " << path << " and " << pe << endl;
					}
				}
			}
		}
	}
	//cout << matrix << endl;
	for (int i = 0; i < paths.size(); i++) {
		//cout << "i: " << i << ", " << paths[i] << endl;
		for (int j = 0; j < paths.size(); j++) {
			//cout << "j: " << j << "/" << paths.size() << ", " << paths[j] << "/" << paths.size() << endl;
			if ((i != j) && (paths[i].canMerge(paths[j]))) {
				//cout << "can merge" << endl;
				bool isLoop     = (paths[i].formLoop(paths[j]));
				PathEnds pe     = paths[i].merge(paths[j]);
				bool insertedPe = false;
				//cout << "new pe " << pe << endl;
				for (const SequencePath &path1: matrix.getPathes(paths[i].getFirst(), paths[i].getLast())) {
					for (const SequencePath &path2: matrix.getPathes(paths[j].getFirst(), paths[j].getLast())) {
						if ((path1.canMerge(path2)) && (! path1.crossesInterior(path2))) {
							SequencePath newPath = path1 + path2;
							//cout << "Adding then " << path1 << "  with " << path2 << " is " << newPath << " and " << pe << endl;
							//cout << "new path " << newPath << endl;
							bool seen = false;
							for (const SequencePath &storedPath: (isLoop)? loops: matrix.getPathes(pe.getFirst(), pe.getLast())) {
								//cout << "stored path " << storedPath << endl;
								if (newPath.contains(storedPath)) {
									//cout << "seen" << endl;
									seen = true;
									break;
								}
							}
							if (! seen) {
								//cout << "adding it" << endl;
								if (isLoop) {
									newPath.setCycle();
									loops.push_back(newPath);
								}
								else {
									matrix.addPath(pe.getFirst(), pe.getLast(), newPath);
									insertedPe = true;
								}
								nbPaths++;
								if (nbPaths % 1000 == 0) {
									cout << "\tFound " << paths.size() << " path ends with " << nbPaths << " paths (now path end " << i << "/" << paths.size() << ").";
									cout << string(70, '\b') << flush;
								}
							}
						}
					}
				}
				if (insertedPe) {
					paths.push_back(pe);
				}
			}
		}
	}
	cout << "\tFound " << paths.size() << " path ends with " << nbPaths << " paths.                                       " << endl;
	for (int i = 0; i < _nodes.size(); i++) {
		for (int j = 0; j < i; j++) {
			for (const SequencePath &path: matrix.getPathes(j, i)) {
				if ((_nodes[path.getNode(0)].isLeaf()) && (_nodes[path.getLastNode()].isLeaf())) {
					_pathes.push_back(path);
				}
			}
		}
	}
	for (SequencePath &loop: loops) {
		_pathes.push_back(loop);
	}
}
*/

/*
void SequenceGraph::findPathes(unsigned int nodeId, mutex &m1) {
	SequencePath path(nodeId);
	SequenceNode &node = _nodes[nodeId];
	node.mark();
	if (isAlone(node)) {
		//cout << "\t\tPath " << node << " is alone" << endl;
		lock_guard<mutex> lock(m1);
		if (checkPath(path)) {
			//cout << "\t\t\tGot it" << endl;
			_pathes.push_back(path);
		}
	}
	else {
		//cout << "\t\tStarting path finding with leaf " << node << " and position " << getOnlyDirection(node) << endl;
		findPathes(path, getOnlyDirection(node), m1);
	}
}

void SequenceGraph::findPathes(SequencePath &path, short position, mutex &m1) {
	int          previousIndex = path.getLastNode();
	SequenceNode previousNode  = _nodes[previousIndex];
	for (int direction = 0; direction < Globals::DIRECTIONS; direction++) {
		short nextPosition = (direction == Globals::DIRECT)? position: 1-position;
		for (int i = 0; i < previousNode.getNbNeighbors(position, direction); i++) {
			unsigned int  nextIndex = previousNode.getNeighbor(position, direction, i);
			SequenceNode &nextNode  = _nodes[nextIndex];
			//cout << "\t\t\tpath is now " << path << ", continuing path finding with node " << nextIndex << " and position " << position << endl;
			if (isLeaf(nextNode)) {
				if (nextIndex < path.getNode(0)) {
					nextNode.mark();
					{
						lock_guard<mutex> lock(m1);
						if (checkPath(path)) {
							_pathes.push_back(path);
						}
					}
				}
			}
			else if (path.contains(nextIndex)) {
				SequencePath cycle(path);
				cycle.trimTo(nextIndex);
				{
					lock_guard<mutex> lock(m1);
					if (checkPath(cycle)) {
						cycle.setCycle();
						_pathes.push_back(cycle);
						//cout << "\t\t\t\tFound cycle " << cycle << endl;
					}
				}
			}
			else {
				nextNode.mark();
				path.addNode(position, direction, nextIndex);
				//cout << "\t\t\t\tGot it!" << endl;
				findPathes(path, nextPosition, m1);
				path.removeLastNode();
			}
		}
	}
}
*/

/*
void SequenceGraph::findCycles() {
	for (SequenceNode &node: _nodes) {
		if (! node.isMarked()) {
			SequencePath path(node.getId());
			findPathes(path, Globals::AFTER);
		}
	}
}

void SequenceGraph::findCycles(SequencePath &path, short position) {
	int          previousIndex = path.getLastNode();
	SequenceNode previousNode  = _nodes[previousIndex];
	for (int direction = 0; direction < Globals::DIRECTIONS; direction++) {
		for (int i = 0; i < previousNode.getNbNeighbors(position, direction); i++) {
			int           nextIndex = previousNode.getNeighbor(position, direction, i);
			SequenceNode &nextNode  = _nodes[nextIndex];
			if (! nextNode.isMarked()) {
				short nextPosition = (direction == Globals::DIRECT)? position: 1-position;
				if (! path.contains(nextIndex)) {
					nextNode.mark();
					path.addNode(nextIndex);
					findCycles(path, nextPosition);
					path.removeLastNode();
				}
				else {
					SequencePath cycle;
					for (int id = path.getSize()-1; path.getNode(id) != nextIndex; id--) {; 
						cycle.addNode(id);
					}
					cycle.setCycle();
					if (checkPath(cycle)) {
						_pathes.push_back(cycle);
					}
				}
			}
		}
	}
}
*/

bool SequenceGraph::checkPath(const SequencePath &path) const {
	for (const SequencePath &storedPath: _pathes) {
		if (path == storedPath) {
			return false;
		}
	}
	return true;
}

void SequenceGraph::selectPathes(const KmerNb threshold) {
	if (_pathes.empty()) {
		return;
	}
	set <unsigned int> nodeIds;
	vector <SequencePath> selectedPathes;
	for (SequencePath &path: _pathes) {
		if (path.getCount() < threshold) {
			path.clear();
		}
	}
	sort(_pathes.begin(), _pathes.end(), PathComparator(*this));
	for (SequencePath &path: _pathes) {
		if (! path.empty()) {
			bool found = false;
			for (unsigned int n = 0; (n < path.getSize()) && (! found); n++) {
				if (nodeIds.find(path.getNode(n)) != nodeIds.end()) {
					found = true;
				}
			}
			if (! found) {
				for (unsigned int n = 0; n < path.getSize(); n++) {
					nodeIds.insert(path.getNode(n));
				}
				selectedPathes.push_back(path);
			}
		}
	}
	_pathes = selectedPathes;
	/*
	cout << "Paths" << endl;
	for (const SequencePath &path: _pathes)
		cout << path << endl;
	*/
}

bool SequenceGraph::comparePathes(const SequencePath &p1, const SequencePath &p2) const {
	return (getCumulatedSize(p1) > getCumulatedSize(p2));
}

vector <SequencePath> &SequenceGraph::getPathes() {
	return _pathes;
}

vector <CountedRepeat> &SequenceGraph::getRepeats() {
	if (_repeats.empty()) {
		for (const SequencePath &path: _pathes) {
			int                  id1       = path.getNode(0);
			const  SequenceNode &n1        = _nodes[id1];
			const  Sequence     &seq1      = n1.getSequence();
			string               str1      = (n1.isDirect())? seq1.getWord(Globals::DIRECT): seq1.getWord(Globals::REVERSE);
			short                position  = Globals::POSITIONS;
			short                direction = Globals::DIRECT;
			//cout << "Making sequence from path " << path << endl;
			//cout << "\tStarting with node " << n1 << " and string " << str1 << endl;
			for (unsigned int j = 1; j < path.getSize(); j++) {
				int                  id2       = path.getNode(j);
				const  SequenceNode &n2        = _nodes[id2];
				short                p         = path.getNodePosition(j);
				short                d         = path.getNodeDirection(j);
				const  Sequence     &seq2      = n2.getSequence();
				if (position == Globals::POSITIONS) {
					position = p;
				}
				if (d == Globals::REVERSE) {
					direction = 1-direction;
				}
				string str2 = seq2.getWord((n2.isDirect())? direction: 1-direction);
				//cout << "\tAdding node " << n2 << " and string " << str2 << endl;
				if (position == Globals::AFTER) {
					str1 += str2.substr(Globals::KMER-1);
				}
				else {
					str1 = str2.substr(0, str2.size() - Globals::KMER + 1) + str1;
				}
				//cout << "\t\tSumming to " << str1 << endl;
			}
			_repeats.push_back(CountedRepeat(str1, path.getCount(), path.isCycle()));
		}
	}
	return _repeats;
}

void SequenceGraph::clear() {
	_nodes.clear();
	_leaves.clear();
	_pathes.clear();
}

bool SequenceGraph::check () const {
	if (Globals::CHECK.size() < Globals::KMER) {
		cout << "\t\tWarning! Cannot check given sequence for its size is less than k-mer size." << endl;
		return false;
	}
	bool returnValue = false;
	for (unsigned int sequencePos = 0; sequencePos < Globals::CHECK.size() - Globals::KMER + 1; sequencePos++) {
		string part  = Globals::CHECK.substr(sequencePos, Globals::KMER);
		for (const SequenceNode &node: _nodes) {
			for (short direction = 0; direction < Globals::DIRECTIONS; direction++) {
				string sequence = node.getSequence().getWord(direction);
				size_t index    = sequence.find(part);
				if (index != string::npos) {
					cout << "\t\t\t\tGot '" << part << "' @ " << index << " in " << node << endl;
					returnValue = true;
				}
			}
		}
	}
	return returnValue;
}

ostream& operator<<(ostream& output, const SequenceGraph& g) {
	output << "Nodes:\n";
	for (const SequenceNode &node: g._nodes) {
		if (node.isSet()) {
			output << node << "\n";
		}
	}
	if (! g._leaves.empty()) {
		output << "Leaves:\n";
		for (unsigned int leafId: g._leaves) {
			const SequenceNode &node = g._nodes[leafId];
			if (node.isSet()) {
				output << node << "\n";
			}
		}
	}
	if (! g._pathes.empty()) {
		output << "Pathes:\n";
		for (const SequencePath &path: g._pathes) {
			output << path << "\n";
		}
	}
	return output;
}
