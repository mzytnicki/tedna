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
#ifndef ASSEMBLER_HPP
#define ASSEMBLER_HPP 1

#include <iostream>
#include "fastxParser.hpp"
#include "simpleKmerCount.hpp"
#include "repeats.hpp"
using namespace std;

class Assembler {

    private:
		int               _insertSize;
		int               _thresholdPc;
		SimpleKmerCount   _kmerCount;
		const char       *_fileName1, *_fileName2, *_outputFileName;
		KmerNb            _threshold;
		Repeats           _repeats;

    public:
        Assembler(const char *fileName1, const char *fileName2, const char *outputFileName, int insertSize, int thresholdPc);
        void assemble();
		void dump();

	private:
		void readFiles();
		void computeDistributions();
		void findRepeats();
		void openLoops();
		void removeInclusions();
		void mergeRepeats();
		void scaffoldRepeats();
		void removeShortRepeats();
		void check(const string &message) const;
};

#endif

