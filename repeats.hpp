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
#ifndef REPEATS_HPP
#define REPEATS_HPP 1

#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include "globals.hpp"
#include "sequence.hpp"
#include "kmer.hpp"
#include "repeatHolder.hpp"
using namespace std;

class CountedRepeat {
	private:
		Sequence _repeat;
		KmerNb   _count;
		bool     _loop;
	public:
		CountedRepeat(const Sequence &repeat, const KmerNb count, const bool loop = false): _repeat(repeat), _count(count), _loop(loop) {}
		const Sequence &getRepeat() const {
			return _repeat;
		}
		const unsigned int getSize() const {
			return _repeat.getSize();
		}
		void setRepeat(const Sequence &repeat) {
			_repeat = repeat;
		}
		void setLoop(const bool loop) {
			_loop = loop;
		}
		bool isLoop() const {
			return _loop;
		}
		void clear () {
			_repeat.clear();
			_count = 0;
		}
		bool empty () const {
			return _repeat.empty();
		}
		KmerNb getCount () const {
			return _count;
		}
		void setCount (const KmerNb count) {
			_count = count;
		}
		string printFasta(int cpt) const {
			stringstream sstm;
			sstm << "te_" << cpt+1 << "_freq:" << _count;
			return _repeat.printFasta(sstm.str());
		}
		bool operator< (const CountedRepeat &cr) const {
			if (_repeat.getUnambiguousSize() > cr._repeat.getUnambiguousSize()) return true;
			if ((_repeat.getUnambiguousSize() == cr._repeat.getUnambiguousSize()) && (_count > cr._count)) return true;
			if ((_repeat.getUnambiguousSize() == cr._repeat.getUnambiguousSize()) && (_count == cr._count) && (_repeat.getFirstWord() < cr._repeat.getFirstWord())) return true;
			return false;
		}
		friend ostream& operator<<(ostream& output, const CountedRepeat& cr) {
			output << cr._repeat << " (" << cr._count << ")";
			return output;
		}

};

typedef vector<CountedRepeat> RepeatsStructure;

class Repeats {

    private:
		RepeatsStructure _repeats;

    public:
        Repeats ();
        void addRepeat (const CountedRepeat &repeat);
        void addRepeat (RepeatHolder &repeatHolder);
        void addRepeat (const Sequence &repeat, KmerNb count, bool loop = false);
        void addRepeat (const string &repeat, KmerNb count, bool loop = false);
        void addRepeats (const Repeats &repeats);
        void removeRepeat (const int i);
		bool isRemoved (const int i) const;
        void changeRepeat (const int i, const Sequence &repeat, KmerNb count);
		void setLoop (const int i, bool loop);
		unsigned int getNbRepeats () const;
		void sort ();
		void cutToThreshold (const KmerNb threshold);
		const CountedRepeat &getRepeat (const int i) const;
		const CountedRepeat &operator[] (const int i) const;
		const bool empty () const;

		void check (const string &message) const;

		void printFasta (ostream& stream = cout);

		friend ostream& operator<<(ostream& output, const Repeats& r);
};

#endif

