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

#include <cmath>
#include <stack>
#include <algorithm>
#include "globals.hpp"
#include "graphTrimmer.hpp"


GraphTrimmer::GraphTrimmer(SequenceGraph &graph): _graph(graph) {}

void GraphTrimmer::trim(bool merge) {
	//cout << "Starting with\n" << _graph;
	bool changed;
	do {
		if (merge) {
			mergeNodes();
		}
		//cout << "Now\n" << _graph;
		pinchBubbles();
		changed = fuseTips();
		if (removeTips()) changed = true;
	}
	while (changed);
	//cout << "Ending with\n" << _graph;
}

/*
void GraphTrimmer::removeHubs() {
	unsigned int dummyIndex = _graph.getSize();
	_graph.addNode(_graph.getSize(), 0);
	SequenceNode &dummyNode = _graph.getNode(dummyIndex);
	dummyNode.unset();
	for (unsigned int i = 0; i < dummyIndex; i++) {
		SequenceNode &node = _graph.getNode(i);
		if (node.isSet()) {
			for (short position = 0; position < Globals::POSITIONS; position++) {
				if (node.isHub(position)) {
					//cout << "Node " << node << ", " << position << " is a hub.";
					for (short direction = 0; direction < Globals::DIRECTIONS; direction++) {
						for (int ni = 0; ni < node.getNbNeighbors(position, direction); ni++) {
							unsigned int neighborIndex = node.getNeighbor(position, direction, ni);
							SequenceNode &neighbor = _graph.getNode(neighborIndex);
							neighbor.updateLinks(i, dummyIndex, true);
						}
					}
					node.removeLinks(position, Globals::DIRECTIONS);
					//cout << "  Now " << node << endl;
				}
			}
		}
	}
	for (unsigned int i = 0; i < dummyIndex; i++) {
		SequenceNode &node = _graph.getNode(i);
		if (node.isSet()) {
			for (short position = 0; position < Globals::POSITIONS; position++) {
				for (short direction = 0; direction < Globals::DIRECTIONS; direction++) {
					bool empty = true;
					for (int ni = 0; (ni < node.getNbNeighbors(position, direction)) && (empty); ni++) {
						if (node.getNeighbor(position, direction, ni) != dummyIndex) {
							empty = false;
						}
					}
					if (empty) {
						node.removeLinks(position, direction);
						//cout << "Updated " << node << endl;
					}
				}
			}
		}
	}
}
*/

void GraphTrimmer::mergeNodes() {
	//cout << "Merging " << _graph.getSize() << " nodes..." << endl;
	//cout << _graph << endl;
	int nbMerges = 0;
	stack <unsigned int> toBeMerged;
	for (unsigned int i = 0; i < _graph.getSize(); i++) {
		if (_graph.getNode(i).isSet()) {
			toBeMerged.push(i);
		}
	}
	while (! toBeMerged.empty()) {
		unsigned int  i  = toBeMerged.top();
		SequenceNode &n1 = _graph.getNode(i);
		toBeMerged.pop();
		if (n1.isSet()) {
			bool merged = false;
			//cout << "Inspecting node " << i << endl;
			for (short position = 0; (position < Globals::POSITIONS) && (! merged); position++) {
				unsigned int j   = -1;
				int   nbNeigbors = 0;
				short direction  = -1;
				for (short d = 0; (d < Globals::DIRECTIONS) && (nbNeigbors < 2); d++) {
					for (int neighbor = 0; (neighbor < n1.getNbNeighbors(position, d)) && (nbNeigbors < 2); neighbor++) {
						int   id2              = n1.getNeighbor(position, d, neighbor);
						const SequenceNode &n2 = _graph.getNode(id2);
						if (n2.isSet()) {
							j = id2;
							direction = d;
							nbNeigbors++;
						}
					}
				}
				if ((nbNeigbors == 1) && (i != j)) {
					SequenceNode &n2 = _graph.getNode(j);
					nbNeigbors = 0;
					short otherPosition = (direction == Globals::DIRECT)? 1-position: position;
					for (short d = 0; (d < Globals::DIRECTIONS) && (nbNeigbors < 2); d++) {
						for (int neighbor = 0; (neighbor < n2.getNbNeighbors(otherPosition, d)) && (nbNeigbors < 2); neighbor++) {
							const SequenceNode &n3 = _graph.getNode(n2.getNeighbor(otherPosition, d, neighbor));
							if (n3.isSet()) {
								nbNeigbors++;
							}
						}
					}
					if (nbNeigbors == 1) {
						bool loop = false;
						for (short d = 0; (d < Globals::DIRECTIONS) && (! loop); d++) {
							for (int neighbor = 0; (neighbor < n2.getNbNeighbors(1-otherPosition, d)) && (! loop); neighbor++) {
								if (n2.getNeighbor(1-otherPosition, d, neighbor) == i) {
									loop = true;
								}
							}
						}
						if (! loop) {
							//cout << "Merging nodes " << i << " and " << j << ": " << n1 << " and " << n2 << " with position " << position << " and direction " << direction;
							n1.merge(position, direction, n2);
							//cout << ", now have " << n1 << endl;
							for (short p = 0; p < Globals::POSITIONS; p++) {
								for (short d = 0; d < Globals::DIRECTIONS; d++) {
									for (int n = 0; n < n2.getNbNeighbors(p, d); n++) {
										unsigned int k = n2.getNeighbor(p, d, n);
										if (k != i) {
											SequenceNode &n3 = _graph.getNode(k);
											if (n3.isSet()) {
												//cout << "\tNodes of " << k << " were : " << _graph.getNode(k);
												n3.updateLinks(j, i, (direction == Globals::DIRECT));
												//cout << " they are now " << _graph.getNode(k) << endl;
											}
										}
									}
								}
							}
							n2.unset();
							nbMerges++;
							toBeMerged.push(i);
							merged = true;
						}
					}
				}
			}
		}
	}
	//cout << "... done (" << nbMerges << " merges)" << endl;
	//cout << _graph << endl;
}

void GraphTrimmer::pinchBubbles() {
	//cout << "Finding bubbles (graph size is " << _graph.getSize() << ")..." << endl;
	if (_graph.getSize() < 3) {
		return;
	}
	int nbBubbles = 0;
	bool bubble;
	for (unsigned int i = 0; i < _graph.getSize(); i++) {
		const SequenceNode &node = _graph.getNode(i);
		do {
			bubble = false;
			if (node.isSet()) {
				short fork = _graph.isVee(node);
				if (fork != 0) {
					short positions = fork-1;
					vector <short> positionArray;
					if (positions == Globals::POSITIONS) {
						positionArray = {Globals::AFTER, Globals::BEFORE};
					}
					else {
						positionArray = {positions};
					}
					SequencePath path(i);
					for (short position: positionArray) {
						//cout << "Starting bubbles with " << i << " " << position << endl;
						if (pinchBubbles(path, position, false)) {
							pinchBubble(path);
							bubble = true;
							nbBubbles++;
						}
					}
				}
			}
		} while (bubble);
	}
	//cout << "... done (" << nbBubbles << " bubbles found)" << endl;
	//cout << _graph << endl;
}

bool GraphTrimmer::pinchBubbles(SequencePath &path, short position, bool change) {
	if (_graph.getCumulatedSize(path) > Globals::BUBBLE_SIZE) {
		return false;
	}
	unsigned int  firstIndex = path.getNode(0);
	unsigned int  lastIndex  = path.getLastNode();
	SequenceNode &lastNode   = _graph.getNode(lastIndex);
	for (short p = 0; p < Globals::DIRECTIONS; p++) {
		if ((p == position) || ((! change) && (path.getSize() > 2))) {
			for (short direction = 0; direction < Globals::DIRECTIONS; direction++) {
				for (int i = 0; i < lastNode.getNbNeighbors(p, direction); i++) {
					unsigned int nextIndex = lastNode.getNeighbor(p, direction, i);
					//cout << "ni: " << nextIndex << ", n-2: " << path.getNode(path.getSize()-2) << endl;
					const SequenceNode &nextNode = _graph.getNode(nextIndex);
					if (nextNode.isSet()) {
						//cout << "\tContinuing bubble with " << path << " and change " << change << ", added " << nextIndex << endl;
						if ((nextIndex == firstIndex) && ((change) || (p != position))) {
							//cout << "\t\tFound bubble in path " << path << endl;
							return true;
						}
						else if (! path.contains(nextIndex)) {
							path.addNode(nextIndex);
							if (pinchBubbles(path, (direction == Globals::DIRECT)? p: 1-p, (change) || (p != position))) {
								return true;
							}
							path.removeLastNode();
						}
					}
				}
			}
		}
	}
	return false;
}

void GraphTrimmer::pinchBubble(SequencePath &path) {
	int    index = 0;
	KmerNb count = -1;
	for (unsigned int i = 0; i < path.getSize(); i++) {
		if (_graph.getNode(path.getNode(i)).getCount() < count) {
			index = i;
		}
	}
	_graph.getNode(path.getNode(index)).unset();
	//cout << "\t\t\tPinching node " << path.getNode(index) << endl;
}

bool GraphTrimmer::fuseTips() {
	bool fused = false;
	for (unsigned int i = 0; i < _graph.getSize(); i++) {
		const SequenceNode &node = _graph.getNode(i);
		if (node.isSet()) {
			int fork =_graph.isVee(node);
			if (fork > 0) {
				short positions = fork - 1;
				vector <short> positionArray;
				if (positions == Globals::POSITIONS) {
					positionArray = {Globals::AFTER, Globals::BEFORE};
				}
				else {
					positionArray = {positions};
				}
				for (short position: positionArray) {
					fused = fused || fuseTips(node, position);
				}
			}
		}
	}
	//cout << "Tip fusing done" << endl;
	//cout << _graph << endl;
	return fused;
}

bool GraphTrimmer::fuseTips(const SequenceNode &node, short position) {
	unsigned int longestSize = 0;
	unsigned int longestNode = 0;
	for (short direction = 0; direction < Globals::DIRECTIONS; direction++) {
		for (int i = 0; i < node.getNbNeighbors(position, direction); i++) {
			int index = node.getNeighbor(position, direction, i);
			const SequenceNode &nextNode = _graph.getNode(index);
			if (nextNode.isSet()) {
				if ((! _graph.isLeaf(nextNode)) || (node == nextNode)) {
					return false;
				}
				if (nextNode.getSize() > longestSize) {
					longestSize = nextNode.getSize();
					longestNode = index;
				}
			}
		}
	}
	//cout << "Fusing tip of " << node << endl;
	bool fused = false;
	for (short direction = 0; direction < Globals::DIRECTIONS; direction++) {
		for (int i = 0; i < node.getNbNeighbors(position, direction); i++) {
			unsigned int index = node.getNeighbor(position, direction, i);
			if (index != longestNode) {
				SequenceNode &nextNode = _graph.getNode(index);
				if (nextNode.isSet()) {
					//cout << "\t\t\tFusing node " << nextNode << endl;
					nextNode.unset();
					fused = true;
				}
			}
		}
	}
	return fused;
}

bool GraphTrimmer::removeTips() {
	bool removed = false;
	for (unsigned int i = 0; i < _graph.getSize(); i++) {
		const SequenceNode &node = _graph.getNode(i);
		if (node.isSet()) {
			int fork =_graph.isFork(node);
			if (fork > 0) {
				short position = fork - 1;
				removed = removed || removeTips(node, position);
			}
		}
	}
	//cout << "Tip removal done" << endl;
	//cout << _graph << endl;
	return removed;
}

bool GraphTrimmer::removeTips(const SequenceNode &node, short position) {
	unsigned int longestSize = 0;
	unsigned int longestNode = 0;
	for (short direction = 0; direction < Globals::DIRECTIONS; direction++) {
		for (int i = 0; i < node.getNbNeighbors(position, direction); i++) {
			int index = node.getNeighbor(position, direction, i);
			const SequenceNode &nextNode = _graph.getNode(index);
			if (nextNode.isSet()) {
				if ((! _graph.isLeaf(nextNode)) || (node == nextNode)) {
					return false;
				}
				if (nextNode.getSize() > longestSize) {
					longestSize = nextNode.getSize();
					longestNode = index;
				}
			}
		}
	}
	for (short direction = 0; direction < Globals::DIRECTIONS; direction++) {
		for (int i = 0; i < node.getNbNeighbors(position, direction); i++) {
			unsigned int index = node.getNeighbor(position, direction, i);
			if ((index != longestNode) && (! _graph.isLeaf(_graph.getNode(index)))) {
				return false;
			}
		}
	}
	//cout << "Removing tip of " << node << endl;
	bool removed = false;
	for (short direction = 0; direction < Globals::DIRECTIONS; direction++) {
		for (int i = 0; i < node.getNbNeighbors(position, direction); i++) {
			unsigned int index = node.getNeighbor(position, direction, i);
			if (index != longestNode) {
				//cout << "\t\t\tRemoving node " << _graph.getNode(index) << endl;
				_graph.getNode(index).unset();
				removed = true;
			}
		}
	}
	return removed;
}

