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
#ifndef GRAPH_TRIMMER_HPP
#define GRAPH_TRIMMER_HPP 1

#include <iostream>
#include <vector>
#include <unordered_map>
#include <set>
#include "globals.hpp"
#include "sequenceGraph.hpp"
using namespace std;


class GraphTrimmer {

    private:
		SequenceGraph &_graph;

    public:
        GraphTrimmer (SequenceGraph &graph);
		void trim (bool merge);
		//void removeHubs ();

	private:
		void mergeNodes ();
		void pinchBubbles ();
		bool pinchBubbles (SequencePath &path, short position, bool change);
		void pinchBubble (SequencePath &path);
		bool fuseTips ();
		bool fuseTips (const SequenceNode &node, short position);
		bool removeTips ();
		bool removeTips (const SequenceNode &node, short position);
};

#endif
