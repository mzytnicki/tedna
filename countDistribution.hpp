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
#ifndef COUNT_DISTRIBUTION_HPP
#define COUNT_DISTRIBUTION_HPP 1

#include <vector>
#include "globals.hpp"
using namespace std;

class CountDistribution {

    protected:
		KmerNb            _maxCount;
		KmerNb            _minCount;
		vector < KmerNb > _countDistribution;

    public:
    CountDistribution ();
		KmerNb getMin () const;
		KmerNb getMax () const;
		void setMin (const KmerNb count);
		void setMax (const KmerNb count);
    void increase(const KmerNb count);
    KmerNb getModeIndex() const;
    KmerNb getThresholdIndex(float threshold) const;
    KmerNb getNbValues() const;
    void clear ();

		friend ostream& operator<<(ostream& output, const CountDistribution& cd);
};

#endif
