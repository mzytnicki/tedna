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
#ifndef KMER_HPP
#define KMER_HPP 1

#include <cstdint>
#include <iostream>
#include "globals.hpp"
#include "sequence.hpp"
#include "kmerCode.hpp"
using namespace std;

class Kmer {

    private:
		Sequence _sequence;
		KmerCode _firstCode, _secondCode;

    public:
		static const KmerCode UNSET;

        Kmer ();
        Kmer (const string &word);
        Kmer (const Sequence &sequence);
        Kmer (const KmerCode code);
        Kmer (const Kmer &k);
		const bool isSet () const;
        const KmerCode getFirstCode  () const;
        const KmerCode getSecondCode () const;
		const Sequence &getSequence () const;
		const string &getFirstWord () const;
		const string &getSecondWord () const;
		const KmerCode getCodeAfter (int code) const;
		const KmerCode getCodeBefore (int code) const;
		const KmerCode getCodeNeighbor (const short code, const short direction) const;
		const KmerCode getCodeNeighbor (int code) const;
		const short compare (const Kmer &k, const short position = Globals::POSITIONS) const;

		friend ostream& operator<<(ostream& output, const Kmer& k);


	private:
		void computeFirstCode ();
		void computeSecondCode ();
		void computeWords ();
};

#endif
