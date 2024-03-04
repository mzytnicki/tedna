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
#ifndef GLOBALS_HPP
#define GLOBALS_HPP

#include <string>
#include <vector>
using namespace std;

typedef unsigned long KmerNb;
typedef unsigned short Penalty;


class Globals {

	public:
		static constexpr int          NB_BLOCKS            = PRE_NB_BLOCKS;
		static constexpr short        NB_BITS_NUCLEOTIDES  = 2;
		static constexpr short        NUCLEOTIDE_MASK      = 3;
		static constexpr short        NB_NUCLEOTIDES       = 4;
		static constexpr short        DIRECT               = 0;
		static constexpr short        REVERSE              = 1;
		static constexpr short        DIRECTIONS           = 2;
		static constexpr short        AFTER                = 0;
		static constexpr short        BEFORE               = 1;
		static constexpr short        POSITIONS            = 2;
		static constexpr uint32_t     KEY_DELETED          = -1;
		static constexpr uint32_t     KEY_EMPTY            = -2;
		static const     string       VERSION;

		static unsigned short KMER;
		static int            NB_THREADS;
		static long           SIZE_THREAD;
		static unsigned long  NB_READS;
		static KmerNb         MIN_COUNT;
		static float          NB_REPETITIONS;
		static float          FREQUENCY_DIFFERENCE;
		static int            MIN_NB_NODES;
		static unsigned int   MAX_NB_NODES;
		static unsigned int   MIN_TE_SIZE;
		static unsigned int   MAX_TE_SIZE;
		static unsigned int   MAX_PATHS;   
		static unsigned int   NB_SMALL_GRAPHS;
		static short          EROSION_STRENGTH;
		static unsigned int   BUBBLE_SIZE;
		static int            MIN_LTR_SIZE;
		static int            MAX_LTR_SIZE;
		static unsigned short SHORT_KMER_SIZE;
		static Penalty        MAX_MERGE_SIZE;
		static Penalty        MIN_MERGE_SIZE;
		static Penalty        PENALTY_SIZE;
		static Penalty        PENALTY_MISMATCH;
		static Penalty        PENALTY_INDEL;
		static Penalty        MAX_PENALTY;
		static float          MIN_IDENTITY;
		static float          MAX_IDENTITY;
		static unsigned int   MERGE_MAX_NB;
		static unsigned int   MERGE_MAX_NODES;
		static unsigned long  MIN_SCAFFOLD_KMERS;
		static unsigned long  MAX_SCAFFOLD_COUNTS;
		static unsigned int   SCAFFOLD_MAX_EV;
		static bool           FASTA_INPUT;
		static string         CHECK;

		static char getComplement(const char c) {
			switch(c) {
				case 'a':
				case 'A':
					return 'T';
				case 'c':
				case 'C':
					return 'G';
				case 'g':
				case 'G':
					return 'C';
				case 't':
				case 'T':
				case 'u':
				case 'U':
					return 'A';
				default:
					return 'N';
			}
		}

		static int getComplementCode(const int i) {
			switch(i) {
				case 0:
					return 3;
				case 1:
					return 2;
				case 2:
					return 1;
				case 3:
					return 0;
				default:
					return 4;
			}
		}

		static string getReverseComplement(const string &s) {
			int    size = s.size();
			string sr(size, 'A');
			for (int i = 0; i < size; i++) {
				sr[i] = getComplement(s[size - 1 - i]);
			}
			return sr;
		}

		static int getCode(const char c) {
			switch(c) {
				case 'a':
				case 'A':
					return 0;
				case 'c':
				case 'C':
					return 1;
				case 'g':
				case 'G':
					return 2;
				case 't':
				case 'T':
				case 'u':
				case 'U':
					return 3;
				default:
					return 4;
			}
		}

		static char getNucleotide(const int i) {
			switch(i) {
				case 0:
					return 'A';
				case 1:
					return 'C';
				case 2:
					return 'G';
				case 3:
					return 'T';
				default:
					return 'N';
			}
		}

		static int localAlignment(const string &firstString, const string &secondString) {
			unsigned int firstSize = firstString.size(), secondSize = secondString.size();
			vector < int > line(firstSize+1);
			int current, tmp, maxValue = 0;
			for (unsigned int j = 1; j <= firstSize; j++) {
				line[j] = 0;
			}
			for (unsigned int i = 1; i <= secondSize; i++) {
				current = 0;
				for (unsigned int j = 1; j <= firstSize; j++) {
					tmp      = line[j];
					line[j]  = max<int>(max<int>(current + ((secondString[i-1]==firstString[j-1])? 1: -1), 0), max<int>(line[j-1]-1, line[j]-1));
					maxValue = max<int>(line[j], maxValue);
					current  = tmp;
				}
			}
			return maxValue;
		}

};

#endif
