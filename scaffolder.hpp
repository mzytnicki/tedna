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
#ifndef SCAFFOLDER_HPP
#define SCAFFOLDER_HPP 1

#include <vector>
#include <array>
#include <map>
#include <tuple>
#include <iostream>
#include <mutex>
#include <thread>
#include "globals.hpp"
#include "repeats.hpp"
#include "fastqParser.hpp"
#include "sequenceGraph.hpp"
using namespace std;



class Scaffolder {

    private:
		map < int, map < int, array < array < vector < int >, Globals::DIRECTIONS >, Globals::POSITIONS > > > _distances;
		//vector < vector < array < array < int, Globals::DIRECTIONS >, Globals::POSITIONS > > > _modes;
		map < string, vector < tuple < int, int, short > > >  _kmers;
		Repeats      &_inputRepeats;
		Repeats       _outputRepeats;
		const char   *_fileName1, *_fileName2;
		unsigned int  _insertSize;
		unsigned int  _maxEvidencesPerNode;

    public:
        Scaffolder (Repeats &r, const char *fileName1, const char *fileName2, unsigned int insertSize);
        void scaffold ();
		Repeats &getRepeats ();

	private:
		void storeKmers ();
		//void buildStructure ();
		void fillStructure ();
		void fillStructure (FastqParser &parser1, FastqParser &parser2, mutex &m);
		//void removeWeakLinks ();
		//void computeMode();
		int computeMode(const vector < int > &distances) const;
		//void checkDistances();
		//tuple <int, int, short, short> findBestScaffold () const;
		void resetMinEvidences ();
		void buildGraphs();
		void findPaths(SequenceGraph &graph, const vector <unsigned int> &translation);
		CountedRepeat mergePath(SequenceGraph &graph, const SequencePath &path, const vector <unsigned int> &translation);
		bool mergeSequences (const int i, const int j, const short position, const short direction, int distance);

		void updateCount (int repeat1, int repeatPos1, int readPos1, short strand1, int lineSize1, int repeat2, int repeatPos2, int readPos2, short strand2, int lineSize2, mutex &m);

		bool checkCell(int i, int j, short position, short direction) const;
		const vector <int>&getCell(int i, int j, short position, short direction);
		void setCell(const int i, const int j, const short position, const short direction, const vector <int>&values);
		void freeCell(const int i, const int j, const short position, const short direction);

		int stitch(const string &seq1, const string &seq2, const vector <int> &v) const;
		pair <int, int> localAlignment(const string &s1, const string &s2, const unsigned int maxStart, const unsigned int maxEnd) const;

		bool aboveEvidenceThreshold (unsigned int threshold) const;

		friend ostream& operator<<(ostream& output, const Scaffolder& s);
};

#endif
