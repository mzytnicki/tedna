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
#ifndef GRAPH_REPEAT_FINDER_HPP
#define GRAPH_REPEAT_FINDER_HPP 1

#include <iostream>
#include "simpleKmerCount.hpp"
#include "repeats.hpp"
#include "sequenceGraph.hpp"
#include "equations.hpp"
using namespace std;

class GraphRepeatFinder {

    private:
		SimpleKmerCount       &_kmerCount;
		KmerNb                 _threshold;
		vector <KmerCode>      _kmers;
		Repeats                _repeats;

    public:
        GraphRepeatFinder (SimpleKmerCount &km, const KmerNb treshold);
        void findRepeats ();
		Repeats &getRepeats ();

	private:
		KmerCode getFirstKmer ();
		void findRepeat (const Kmer &findKmer);
		void fillFirstGraph (SequenceGraph &graph, const Kmer &firstKmer);
		//SequenceGraph findBestGraph (SequenceGraph &firstGraph);
		void build(SequenceGraph &graph);
		void removeKmers ();
		//void decreaseKmers (Graph &graph);
		void solveProblem(SequenceGraph &graph);
		void addRepeats(SequenceGraph &graph);
		//void gatherRepeats();
};

#endif
