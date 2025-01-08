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

#include <ostream>
#include <iostream>
#include "globals.hpp"
#include "countDistribution.hpp"
using namespace std;

CountDistribution::CountDistribution (): _maxCount(0), _minCount(0) { }

KmerNb CountDistribution::getMin() const {
  return _minCount;
}

KmerNb CountDistribution::getMax() const {
  return _maxCount;
}

void CountDistribution::setMin(const KmerNb count) {
  _minCount = count;
}

void CountDistribution::setMax(const KmerNb count) {
  _maxCount = count;
	_countDistribution.assign(count+1, 0);
}

void CountDistribution::increase(const KmerNb count) {
	++_countDistribution[count];
}

KmerNb CountDistribution::getModeIndex() const {
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

KmerNb CountDistribution::getThresholdIndex(float threshold) const {
  KmerNb sum = 0;
  for (KmerNb i = _maxCount; i != 0; i--) {
    sum += i * _countDistribution[i];
    if (sum >= threshold) {
      return i;
    }
  }
  return 0;
}

void CountDistribution::clear() {
	_countDistribution.clear();
}

ostream& operator<<(ostream& output, const CountDistribution& cd) {
	for (KmerNb nb = 0; nb <= cd._maxCount; ++nb) {
		if (cd._countDistribution[nb] != 0) {
			cout << "\t\t" << nb << ": " << cd._countDistribution[nb] << endl;
		}
	}
	return output;
}
