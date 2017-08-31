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

#include <limits>
#include <vector>
#include "globals.hpp"
#include "fastqParser.hpp"

FastqParser::FastqParser(const char *fileName): _file(fileName), _pos(-1), _over(false), _allRead(false), _end(numeric_limits<unsigned long long>::max()), _lineNb(0), _readId(0) {
	if (! _file.is_open()) {
		throw "Error! Input file '" + string(fileName) + "' cannot be opened!";
	}
}

unsigned long FastqParser::getReadId() const {
	return _readId;
}

void FastqParser::reset() {
	_readId       = 0;
	_line         = "";
	_over         = false;
	_allRead      = false;
	_pos          = -1;
	_word         = "";
	_end          = numeric_limits<unsigned long long>::max();
	_file.clear();
	_file.seekg(0, ios::beg);
	_sequence.clear();
}

void FastqParser::goTo(unsigned long long start, bool skip) {
	reset();
	_file.seekg(start, ios::beg);
	if (! _file.good()) {
		//cout << "file is not good after go to" << endl;
		_allRead = true;
		_over    = true;
		return;
	}
	if (skip) {
		if (! goToNextLine()) {
			//cout << "file is not good after next line" << endl;
			_allRead = true;
			_over    = true;
		}
	}
}

void FastqParser::endTo(unsigned long long end) {
	_end = end;
	//cout << this << ": end is " << _end << endl;
}

bool FastqParser::isOver() const {
	return _over;
}

bool FastqParser::isAllRead() const {
	return _allRead;
}

void FastqParser::getNextKmer() {
	findNextCheckedKmer();
}

void FastqParser::getNextLine() {
	readNewLine();
}

string &FastqParser::getLine() {
	if (_line.empty()) {
		getNextLine();
	}
	return _line;
}

string &FastqParser::getWord() {
	if (_word.empty()) {
		getNextKmer();
	}
	return _word;
}

Sequence &FastqParser::getSequence() {
	if (_sequence.empty()) {
		getNextKmer();
	}
	return _sequence;
}

void FastqParser::findNextKmer() {
	if (_word.length() == Globals::KMER) {
		_word = _word.substr(1, Globals::KMER-1);
	}
	while (_word.length() < Globals::KMER) {
		if ((_pos >= _line.length()) || (_pos == static_cast<unsigned int>(-1))) {
			readNewLine();
			//cout << this << ": pos is " << _file.tellg() << ", end is " << _end << endl;
			if (isOver()) {
				//cout << "normal end" << endl;
				return;
			}
			_word = _line.substr(0, Globals::KMER);
			_pos  = Globals::KMER;
		}
		else {
			_word += _line[_pos];
			_pos++;
		}
		_sequence = Sequence(_word);
	}
}

void FastqParser::findNextCheckedKmer() {
	do {
		findNextKmer();
	}
	//while ((not isOver()) && (_sequence.isAmbiguous()));
	while ((not isOver()) && ((_sequence.isAmbiguous()) || (_sequence.isLowComplexity())));
}

void FastqParser::readNewLine() {
	do {
		getline(_file, _line);
		if (! _file.good()) {
			//cout << "read new line end" << endl;
			_allRead = true;
			_over    = true;
			return;
		}
		_lineNb++;
	}
	while (_lineNb % 4 != 2);
	_pos  = 0;
	_word = "";
	_readId++;
	if (_file.tellg() > _end) {
		//cout << "end of my part" << endl;
		_over = true;
	}
}

bool FastqParser::goToNextLine() {
	vector<string> lines(4);
	unsigned long long pos = _file.tellg();
	for (int i = 0; i < 4; i++) {
		getline(_file, lines[i]);
		if (i > 0) {
			if (! _file.good()) {
				//cout << "file is not good after go to fastq" << endl;
				return false;
			}
		}
	}
	int start;
	if ((lines[1][0] == '@') && (lines[3][0] == '+')) {
		start = 1;
	}
	else if (lines[2][0] == '+') {
		start = 4;
	}
	else if (lines[1][0] == '+') {
		start = 3;
	}
	else {
		start = 2;
	}
	_file.seekg(pos);
	_lineNb = 0;
	//cout << "start is " << start << endl;
	for (int i = 0; i < start; i++) {
		getline(_file, lines[0]);
		//cout << lines[0] << endl;
	}
	return true;
}
