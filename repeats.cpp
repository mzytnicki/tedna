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

#include <algorithm>
#include "repeats.hpp"

Repeats::Repeats() { }


void Repeats::addRepeat(const CountedRepeat &repeat) {
	_repeats.push_back(repeat);
}

void Repeats::addRepeat(RepeatHolder &repeatHolder) {
	addRepeat(CountedRepeat(repeatHolder.getRepeat(), repeatHolder.getAverage()));
}

void Repeats::addRepeat(const Sequence &repeat, KmerNb count, bool loop) {
	addRepeat(CountedRepeat(repeat, count, loop));
}

void Repeats::addRepeat(const string &repeat, KmerNb count, bool loop) {
	Sequence sequence(repeat);
	addRepeat(sequence, count, loop);
}

void Repeats::addRepeats(const Repeats &repeats) {
	for (const CountedRepeat &cr: repeats._repeats) {
		addRepeat(cr);
	}
}

void Repeats::changeRepeat(const int i, const Sequence &repeat, KmerNb count) {
	_repeats[i] = CountedRepeat(repeat, count);
}

void Repeats::setLoop(const int i, bool loop) {
	_repeats[i].setLoop(loop);
}

void Repeats::removeRepeat(const int i) {
	_repeats[i].clear();
}

bool Repeats::isRemoved(const int i) const {
	return _repeats[i].empty();
}

unsigned int Repeats::getNbRepeats() const {
	return _repeats.size();
}

void Repeats::sort() {
	RepeatsStructure repeats;
	for (const CountedRepeat &repeat: _repeats) {
		if (repeat.getCount() > 0) {
			repeats.push_back(repeat);
		}
	}
	std::sort(repeats.begin(), repeats.end());
	_repeats = repeats;
}

void Repeats::cutToThreshold(const KmerNb threshold) {
	RepeatsStructure repeats;
	for (CountedRepeat &cr: _repeats) {
		if (cr.getCount() >= threshold) {
			repeats.push_back(cr);
		}
	}
	_repeats = repeats;
}

const CountedRepeat &Repeats::getRepeat(const int i) const {
	return _repeats[i];
}

const CountedRepeat &Repeats::operator[] (const int i) const {
	return getRepeat(i);
}

const bool Repeats::empty () const {
	return _repeats.empty();
}

void Repeats::printFasta (ostream& stream) {
	sort();
	for (unsigned int i = 0; i < _repeats.size(); i++) {
		if (! _repeats[i].empty()) {
			stream << _repeats[i].printFasta(i);
		}
	}
}

void Repeats::check (const string &message) const {
	if (Globals::CHECK.empty()) {
		return;
	}
	cout << "\t\t" << message << endl;
	for (const CountedRepeat &repeat: _repeats) {
		for (short direction = 0; direction < Globals::DIRECTIONS; direction++) {
			string sequence = repeat.getRepeat().getWord(direction);
			int score = Globals::localAlignment(Globals::CHECK, sequence);
			if (static_cast<unsigned int>(score) >= Globals::CHECK.size() / 2) {
				cout << "\t\t\tScore of " << score << " with " << sequence << endl;
			}
		}
	}
}

ostream& operator<<(ostream& output, const Repeats& r) {
	for (unsigned int i = 0; i < r._repeats.size(); i++) {
		output << i << ": " << r._repeats[i] << endl;
	}
	return output;
}

