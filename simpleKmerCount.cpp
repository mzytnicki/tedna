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
#include "globals.hpp"
#include "simpleKmerCount.hpp"

static const KmerNb MAX_KMER_NB = numeric_limits<KmerNb>::max();

SimpleKmerCount::SimpleKmerCount(): _maxCount(0), _minCount(0), _nbValues(0) { }

void SimpleKmerCount::addKmer(const Kmer &kmer, bool insert) {
	//cout << "Adding " << kmer.getFirstCode() << endl;
	++_map[kmer.getFirstCode()];
}

void SimpleKmerCount::addKmer(const Kmer &kmer, mutex &m) {
	lock_guard<mutex> lock(m);
	addKmer(kmer);
}

KmerNb SimpleKmerCount::getCount(const Kmer &kmer) const {
	auto it = _map.find(kmer.getFirstCode());
	if (it == _map.end()) {
		return 0;
	}
	return it->second;
}

bool SimpleKmerCount::isPresent(const KmerCode &kmerCode) const {
	return (_map.find(kmerCode) != _map.end());
}

bool SimpleKmerCount::isPresent(const Kmer &kmer) const {
	return (isPresent(kmer.getFirstCode()));
}

void SimpleKmerCount::decreaseNb(const KmerCode &kmerCode, const KmerNb nb) {
	if (! isPresent(kmerCode)) {
		return;
	}
	KmerNb count = _map[kmerCode];
	if (count <= nb + _minCount) {
		_nbValues -= count;
		_map.erase(kmerCode);
	}
	else {
		_nbValues -= nb;
		_map[kmerCode] -= nb;
	}
}

void SimpleKmerCount::decreaseNb(const Kmer &kmer, const KmerNb nb) {
	return decreaseNb(kmer.getFirstCode(), nb);
}

void SimpleKmerCount::remove(const KmerCode &kmerCode) {
	_nbValues -= _map[kmerCode];
	_map.erase(kmerCode);
}

void SimpleKmerCount::remove(const Kmer &kmer) {
	KmerCode code = kmer.getFirstCode();
	_nbValues -= _map[code];
	_map.erase(code);
}

void SimpleKmerCount::computeCountDistribution() {
	_maxCount = 0;
	_nbValues = 0;
	for (auto it = _map.begin(); it != _map.end(); ++it) {
		_maxCount = max<KmerNb>(_maxCount, it->second);
		_nbValues += it->second;
	}
	_countDistribution.assign(_maxCount+1, 0);
	for (auto it = _map.begin(); it != _map.end(); ++it) {
		++_countDistribution[it->second];
	}
}

void SimpleKmerCount::printCountDistribution() const {
	for (KmerNb nb = 0; nb <= _maxCount; ++nb) {
		if (_countDistribution[nb] != 0) {
			cout << "\t\t" << nb << ": " << _countDistribution[nb] << endl;
		}
	}
}

void SimpleKmerCount::setMinCount(const KmerNb count) {
	_minCount = count;
}

KmerNb SimpleKmerCount::getMaxCount() const {
	return _maxCount;
}

KmerNb SimpleKmerCount::getMaxCountDistribution() const {
	KmerNb index = 0;
	KmerNb value = 0;
	for (KmerNb i = Globals::MIN_COUNT; i <= _maxCount; i++) {
		if (_countDistribution[i] > value) {
			index = i;
			value = _countDistribution[i];
		}
	}
	return index;
}

KmerNb SimpleKmerCount::getThreshold(int percent) const {
	KmerNb sum = 0;
	KmerNb threshold = _nbValues * percent / 100;
	for (KmerNb i = _maxCount; i != 0; i--) {
		sum += i * _countDistribution[i];
		if (sum >= threshold) {
			return i;
		}
	}
	return 0;
}

void SimpleKmerCount::removeUnder(KmerNb nb) {
	for (auto it = _map.begin(); it != _map.end(); ) {
		if (it->second < nb) {
			_nbValues -= it->second;
			_map.erase(it++);
		}
		else {
			++it;
		}
	}
}

KmerCode SimpleKmerCount::getMostFrequent() {
	KmerCode index = 0;
	KmerNb   value = 0;
	for (auto it = _map.begin(); it != _map.end(); ++it) {
		if (it->second > value) {
			index = it->first;
			value = it->second;
		}
	}
	return index;
}

KmerCode SimpleKmerCount::getLeastFrequent() {
	KmerCode index = 0;
	KmerNb   value = -1;
	for (auto it = _map.begin(); it != _map.end(); ++it) {
		if (it->second < value) {
			index = it->first;
			value = it->second;
		}
	}
	return index;
}

pair <KmerCode, KmerNb> SimpleKmerCount::getRandom() {
	auto it = _map.begin();
	if (it == _map.end()) {
		return make_pair(Kmer::UNSET, 0);
	}
	return *it;
}

bool SimpleKmerCount::empty() const {
	return _map.empty();
}

void SimpleKmerCount::clear() {
	_map.clear();
	_countDistribution.clear();
}

unsigned int SimpleKmerCount::getSize() const {
	return _map.size();
}

ostream& operator<<(ostream& output, SimpleKmerCount& kc) {
	for (auto it = kc._map.begin(); it != kc._map.end(); ++it) {
		output << Kmer(it->first) << "\t" << it->second << endl;
	}
	return output;
}
