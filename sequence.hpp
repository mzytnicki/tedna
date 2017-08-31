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
#ifndef SEQUENCE_HPP
#define SEQUENCE_HPP 1

#include <string>
#include <iostream>
#include <tuple>
#include "globals.hpp"
using namespace std;

class Sequence {

    private:
		string _firstWord;
		string _secondWord;

    public:
        Sequence ();
        Sequence (const string &word);
        Sequence (const Sequence &s);
		const bool empty() const;
		const unsigned int getSize () const;
		const unsigned int getUnambiguousSize () const;
		const string &getFirstWord ()  const;
		const string &getSecondWord () const;
		const string &getWord (const short i) const;
		void setFirstWord (const string &word);
		bool addFront (const string &nucleotides);
		bool addBack  (const string &nucleotides);
		tuple<short, short, bool> merge (const Sequence &k, short direction = Globals::DIRECTIONS, int size = -1);
		bool addCode (const int code);
		void clear ();
		bool updateWords();
		const short compareNext (const Sequence &s, const short position = Globals::POSITIONS) const;

		bool isAmbiguous () const;
		bool isLowComplexity () const;

		string printFasta(string name) const;

		friend bool operator==(const Sequence &s1, const Sequence &s2);
		friend ostream& operator<<(ostream& output, const Sequence& s);


	private:
		bool computeSecondWord ();
};

#endif
