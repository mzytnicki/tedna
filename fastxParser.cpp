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
#include "fastxParser.hpp"

FastxParser::FastxParser(unsigned int b, unsigned int s, const char *fileName): _file(fileName), _pos(-1), _over(false), _allRead(false), _end(numeric_limits<unsigned long long>::max()), _lineNb(0), _readId(0), _blockSize(b), _sequenceLine(s) {
	if (! _file.is_open()) {
		throw "Error! Input file '" + string(fileName) + "' cannot be opened!";
	}
}

unsigned long FastxParser::getReadId() const {
	return _readId;
}

void FastxParser::reset() {
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

void FastxParser::goTo(unsigned long long start, bool skip) {
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

void FastxParser::endTo(unsigned long long end) {
	_end = end;
	//cout << this << ": end is " << _end << endl;
}

bool FastxParser::isOver() const {
	return _over;
}

bool FastxParser::isAllRead() const {
	return _allRead;
}

void FastxParser::getNextKmer() {
	findNextCheckedKmer();
}

void FastxParser::getNextLine() {
	readNewLine();
}

string &FastxParser::getLine() {
	if (_line.empty()) {
		getNextLine();
	}
	return _line;
}

string &FastxParser::getWord() {
	if (_word.empty()) {
		getNextKmer();
	}
	return _word;
}

Sequence &FastxParser::getSequence() {
	if (_sequence.empty()) {
		getNextKmer();
	}
	return _sequence;
}

void FastxParser::findNextKmer() {
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

void FastxParser::findNextCheckedKmer() {
	do {
		findNextKmer();
	}
	//while ((not isOver()) && (_sequence.isAmbiguous()));
	while ((not isOver()) && ((_sequence.isAmbiguous()) || (_sequence.isLowComplexity())));
}

void FastxParser::readNewLine() {
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
	while (_lineNb % _blockSize != _sequenceLine);
	_pos  = 0;
	_word = "";
	_readId++;
	if (_file.tellg() > _end) {
		//cout << "end of my part" << endl;
		_over = true;
	}
}

FastqParser::FastqParser (const char *fileName): FastxParser(4, 2, fileName) { }

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

FastaParser::FastaParser (const char *fileName): FastxParser(2, 0, fileName) { }

bool FastaParser::goToNextLine() {
	string line;
	unsigned long long pos = _file.tellg();
	getline(_file, line);
	if (! _file.good()) {
		//cout << "file is not good after go to fastq" << endl;
		return false;
	}
	//cout << line << endl;
	if (line[0] == '>') {
		_file.seekg(pos);
	}
	return true;
}
