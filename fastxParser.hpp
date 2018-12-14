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
#ifndef FASTX_PARSER_HPP
#define FASTX_PARSER_HPP 1

#include <string>
#include <iostream>
#include <fstream>
#include "sequence.hpp"
using namespace std;

class FastxParser {

	protected:
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
		unsigned int       _blockSize;
		unsigned int       _sequenceLine;

	public:
		FastxParser (unsigned int b, unsigned int s, const char *fileName);
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
		virtual ~FastxParser() {}

	protected:
		void findNextKmer ();
		void findNextCheckedKmer ();
		void readNewLine ();
		virtual bool goToNextLine () = 0;
};

class FastqParser: public FastxParser  {

	public:
		FastqParser (const char *fileName);

	protected:
		virtual bool goToNextLine ();
};

class FastaParser: public FastxParser  {

	public:
		FastaParser (const char *fileName);

	protected:
		virtual bool goToNextLine ();
};

#endif


