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
#ifndef FASTQ_PARSER_HPP
#define FASTQ_PARSER_HPP 1

#include <string>
#include <iostream>
#include <fstream>
#include "sequence.hpp"
using namespace std;

class FastqParser {

	private:
		ifstream           _file;
		unsigned int       _pos;
		bool               _over;
		bool               _allRead;
		long long          _end;
		unsigned long      _lineNb;
		unsigned long      _readId;
		string             _line;
		string             _word;
		Sequence           _sequence;

    public:
        FastqParser (const char *fileName);
		unsigned long getReadId() const;
		bool isOver () const;
		bool isAllRead () const;
		void getNextKmer ();
		void getNextLine ();
		string &getLine ();
		string &getWord ();
		Sequence &getSequence ();
		void reset ();
		void goTo (unsigned long long start, bool skip = false);
		void endTo (unsigned long long end);

	protected:
		void findNextKmer ();
		void findNextCheckedKmer ();
		void readNewLine ();
		bool goToNextLine ();
};

#endif


