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
#ifndef LOOP_OPENER_HPP
#define LOOP_OPENER_HPP 1

#include <iostream>
#include <vector>
#include <map>
#include <mutex>
#include "globals.hpp"
#include "kmerCode.hpp"
#include "kmer.hpp"
#include "sequence.hpp"
#include "repeats.hpp"
using namespace std;


class LoopOpener {

    private:
		const char *_fileName;
		Repeats     _repeats;
		int         _nbLoops;

		map < KmerCode, KmerNb > _count;

    public:
        LoopOpener (const char *fileName, const Repeats &repeats);
		void openLoops ();
		const Repeats &getRepeats();

		friend ostream& operator<<(ostream& output, const LoopOpener& lo);

	private:
		void reset ();
		void addLoop (const Sequence &loop);
		void readReads ();

		pair <int, int> getOverRepresented (const string &loop);
		int getWeakPosition (const string &loop);

		bool empty() const;
};

#endif

