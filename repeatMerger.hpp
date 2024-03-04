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
#ifndef REPEAT_MERGER_HPP
#define REPEAT_MERGER_HPP 1

#include <array>
#include <vector>
#include <iostream>
#include "globals.hpp"
#include "repeats.hpp"
#include "kmerSet.hpp"
#include "sequenceComparator.hpp"
#include "sequenceGraph.hpp"
using namespace std;


class ComparisonData {

	private:
		bool    _set;
		Penalty _score;
		float   _identity;
		int     _startFirst, _endSecond;

	public:
		ComparisonData(): _set(false) {}

		ComparisonData(const SequenceComparator &sc, string &first, string &second):
			_set(true),
			_score(sc.getBestScore()),
			_identity(sc.getIdentity()),
			_startFirst(sc.getStartFirst()),
			_endSecond(sc.getEndSecond()) {}

		void unset () {
			_set = false;
		}

		bool isSet () const {
			return _set;
		}

		Penalty getScore () const {
			return _score;
		}

		float getIdentity () const {
			return _identity;
		}

		int getStartFirst () const {
			return _startFirst;
		}

		int getEndSecond () const {
			return _endSecond;
		}

		void copy (const ComparisonData &cd, const short direction, const short position) {
			//cout << "    A " << *this << "  " << cd << endl;
			_set = cd._set;
			if (! _set) {
				return;
			}
			_score    = cd._score;
			_identity = cd._identity;
			//cout << "    B " << *this << endl;
		}

		friend ostream& operator<<(ostream& output, const ComparisonData& cd) {
			if (! cd.isSet()) {
				output << "(unset)";
			}
			else {
				output << cd._score << ", " << cd._identity << ": (" << cd._startFirst << ", " << cd._endSecond << ")";
			}
			return output;
		}
};


class RepeatMerger {

	private:
		unsigned int _size;
		KmerNb       _minCount;
		KmerSet*     _kmerSets;
		Repeats      _inputRepeats;
		Repeats      _outputRepeats;
		Penalty      _maxPenalty;
		vector < vector < array < array < ComparisonData, Globals::POSITIONS > , Globals::DIRECTIONS > > > _comparisons;

	public:
		RepeatMerger (const Repeats &repeats, const KmerNb minCount = 0);
		void addRepeat (const Sequence &sequence, const KmerNb count);
		void mergeRepeats ();
		void report ();
		Repeats &getRepeats();

	private:
		void buildStructure ();
		void cleanStructure ();
		void fillStructure ();
		void buildGraphs();
		void findPaths(SequenceGraph &graph, const vector <unsigned int> &translation);
		//void findBestMerge ();
		//void mergeSequences ();
		CountedRepeat mergePath (SequenceGraph &graph, const SequencePath &path, const vector <unsigned int> &translation);

		void resetMaxPenalty ();
		bool underPenaltyThreshold (unsigned int threshold) const;

		const ComparisonData &getComparison(int i, int j, short direction, short position) const;
		bool isSet(unsigned int i, unsigned int j, short k, short l) const;
		void setCell(unsigned int i, unsigned int j, short k, short l, SequenceComparator &sc, string &first, string &second);

		void compare(SequenceComparator &comparator, unsigned int i, unsigned int j);

		friend ostream& operator<<(ostream& output, const RepeatMerger& rm);
};

#endif
