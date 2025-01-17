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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string>
#include <iostream>
#include "globals.hpp"
#include "kmerParser.hpp"
using namespace std;

KmerParser::KmerParser (const char *fileName): _file(fileName), _maxCount(0), _nbValues(0) { }

KmerNb KmerParser::getMaxCount () const {
  return _maxCount;
}

KmerNb KmerParser::getNbValues () const {
  return _nbValues;
}

void KmerParser::getStats () {
  _file.clear();
  _file.seekg(0);
  for (string line; getline(_file, line); ++_nbValues) {
    if (! line.empty()) {
      if (line[0] == '>') {
        KmerNb nbInt = atoi(&line[1]);
        _maxCount = max(_maxCount, nbInt);
      }
    }
  }
}

void KmerParser::fill (CountDistribution &cd) {
  _file.clear();
  _file.seekg(0);
  KmerNb cpt = 0;
  cout << "Reading k-mer file...\n";
  for (string line; getline(_file, line); ++cpt) {
    if (! line.empty()) {
      if (line[0] == '>') {
        KmerNb nbInt = atoi(&line[1]);
        cd.increase(nbInt);
      }
    }
    if (cpt % 100000000 == 0) cout << "\t" << cpt << " lines read\n";
  }
  cout << "\t" << cpt << " lines read, done.\n";
}

/*
void KmerParser::fillKmers (BBhashStr &bb, KmerNb minCount) {
  _file.clear();
  _file.seekg(0);
  KmerNb nb = 0;
  KmerNb cpt = 0;
  KmerNb cptAdded = 0;
  cout << "Filling k-mer table...\n";
  for (string line; getline(_file, line); ++cpt) {
    if (! line.empty()) {
      if (line[0] == '>') {
        nb = atoi(&line[1]);
      }
      else {
        if (nb >= minCount) {
          Kmer k(line);
          bb.add(k.getFirstCode());
          ++cptAdded;
        }
      }
    }
    if (cpt % 100000000 == 0) cout << "\t" << cpt << " lines read\n";
  }
  cout << "\t" << cpt << " lines read, " << cptAdded << " k-mers inserted, done.\n";
}

void KmerParser::fillCounts (BBhashStr &bb, KmerNb minCount) {
  _file.clear();
  _file.seekg(0);
  KmerNb nb = 0;
  KmerNb cpt = 0;
  cout << "Adding k-mer counts...\n";
  for (string line; getline(_file, line); ++cpt) {
    if (! line.empty()) {
      if (line[0] == '>') {
        nb = atoi(&line[1]);
      }
      else {
        if (nb >= minCount) {
          Kmer k(line);
          bb.setCount(k.getFirstCode(), nb);
        }
      }
    }
    if (cpt % 100000000 == 0) cout << "\t" << cpt << " lines read\n";
  }
  cout << "\t" << cpt << " lines read, done.\n";
}
*/

void KmerParser::fill(hash_t &map, KmerNb minCount) {
  _file.clear();
  _file.seekg(0);
  KmerNb nb = 0;
  KmerNb cpt = 0;
  cout << "Adding k-mer counts...\n";
  for (string line; getline(_file, line); ++cpt) {
    if (! line.empty()) {
      if (line[0] == '>') {
        nb = atoi(&line[1]);
      }
      else {
        if (nb >= minCount) {
          Kmer k(line);
          map[k.getFirstCode()] = nb;
        }
      }
    }
    if (cpt % 100000000 == 0) cout << "\t" << cpt << " lines read\n";
  }
  cout << "\t" << cpt << " lines read, done.\n";
}
