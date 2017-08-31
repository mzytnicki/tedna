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
#ifndef SEQUENCE_GRAPH_HPP
#define SEQUENCE_GRAPH_HPP

#include <iostream>
#include <vector>
#include <set>
#include <algorithm>
#include <mutex>
#include "globals.hpp"
#include "sequence.hpp"
#include "repeats.hpp"
using namespace std;


class SequenceNode {

	private:
		unsigned int _id;
		KmerNb       _count;
		Sequence     _sequence;
		bool         _set;
		bool         _marked;
		bool         _direct;
		vector <unsigned int> _neighbors[Globals::POSITIONS][Globals::DIRECTIONS];

	public:
		SequenceNode (unsigned int id = -1, KmerNb count = 0);
		SequenceNode (unsigned int id, KmerNb count, const Sequence &sequence);

		void addNeighbor(unsigned int id, short position, short direction);

		void setId(unsigned int id);
		unsigned int getId() const;
		unsigned int getSize() const;
		const Sequence &getSequence() const;
		void setSequence(const Sequence sequence);
		KmerNb getCount() const;
		void unset();
		bool isSet() const;
		void mark();
		bool isMarked() const;
		bool isDirect() const;
		bool isLeaf() const;
		short getOnlyDirection() const;
		int getNbNeighbors(short position, short direction) const;
		const unsigned int getNeighbor(short position, short direction, int i) const;

		bool merge (short position, short direction, const SequenceNode &n);
		void updateLinks(unsigned int oldId, unsigned int newId, bool direct);
		void removeLinks(short position, short direction);

		friend bool operator==(const SequenceNode& n1, const SequenceNode& n2);
		friend ostream& operator<<(ostream& output, const SequenceNode& n);
};


class SequencePath {

	private:
		vector < tuple <short, short, unsigned int> > _nodes;
		set    < unsigned int > _nodeIds;
		bool   _cycle;
		KmerNb _count;

	public:
		SequencePath ();
		SequencePath (const SequencePath &path);
		SequencePath (unsigned int id);
		SequencePath (unsigned int id1, short position, short direction, unsigned int id2);

		void addNode(unsigned int id);
		void addNode(short position, short direction, unsigned int id);
		void removeLastNode();

		unsigned int getSize() const;
		unsigned int getNode(int i) const;
		int getNodePosition(int i) const;
		int getNodeDirection(int i) const;
		unsigned int getLastNode() const;

		bool empty() const;
		void clear();
		void reverse();
		
		void setCycle();
		bool isCycle() const;

		void setCount(const KmerNb count);
		KmerNb getCount() const;

		bool contains(int id) const;
		bool contains(const SequencePath &p) const;
		bool crosses(const SequencePath &p) const;
		bool crossesInterior(const SequencePath &p) const;
		bool canMerge(const SequencePath &p) const;
		bool formCycle(const SequencePath &p) const;

		void trimTo(const unsigned int destination);

		vector < SequencePath > addNode(const SequenceNode &linkNode, const SequenceNode &newNode) const;
		SequencePath merge(const SequenceNode &linkNode, const SequencePath &path) const;

		friend SequencePath operator+(const SequencePath& p1, const SequencePath& p2);
		friend bool operator==(const SequencePath& p1, const SequencePath& p2);
		friend ostream& operator<<(ostream& output, const SequencePath& p);
};


class SequenceGraph {

	private:
		unsigned int _maxPaths;
		vector <SequenceNode> _nodes;
		vector <unsigned int>  _leaves;
		vector <SequencePath>  _pathes;
		vector <CountedRepeat> _repeats;

	public:
		SequenceGraph ();
		SequenceGraph (int size);

		void setMaxPaths(const unsigned int maxPaths);

		unsigned int getSize() const;
		bool isSmall() const;
		bool isBig() const;

		void addNode(unsigned int id, KmerNb count);
		void addNode(unsigned int id, KmerNb count, const Sequence &sequence);
		void addLink(unsigned int n1, short pos1, short dir1, unsigned int n2);
		SequenceNode &getNode(const unsigned int i);
		unsigned int getCumulatedSize(const SequencePath &path) const;
		bool findAllPathes();
		void findGreedyPathes();
		void selectPathes(const KmerNb count = 0);
		vector <SequencePath> &getPathes();
		vector <CountedRepeat> &getRepeats();

		void clear();

		bool isLeaf(const SequenceNode &n) const;
		int  isVee(const SequenceNode &n) const;
		int  isFork(const SequenceNode &n) const;
		bool isAlone(const SequenceNode &n) const;

		bool check () const;

		friend ostream& operator<<(ostream& output, const SequenceGraph& g);

	private:
		void findLeaves();
		bool findPathes();
		//void findCycles();
		//void findPathes(unsigned int leafId, mutex &m1);
		//void findPathes(SequencePath &path, short position, mutex &m1);
		bool checkPath(const SequencePath &p) const;
		bool comparePathes(const SequencePath &p1, const SequencePath &p2) const;
		void countNeighbors(const SequenceNode &node, KmerNb *count) const;
		short getOnlyDirection(const SequenceNode &n) const;
		void removeMarkedNodes();

		unsigned int findMostSeenNode() const;
		pair < unsigned int, bool > findMostSeenNeighbor(const unsigned int nodeId, const unsigned short position) const;

		struct PathComparator { 
			const SequenceGraph& _graph;
			PathComparator(const SequenceGraph& graph) : _graph(graph) {}
			bool operator()(const SequencePath &p1, const SequencePath &p2) { 
				return (_graph.comparePathes(p1, p2));
			}
		};
};

#endif
