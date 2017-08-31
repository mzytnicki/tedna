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
#ifndef EQUATIONS_HPP
#define EQUATIONS_HPP 1

#include <cmath>
#include <iostream>
#include <map>
#include <vector>
using namespace std;


class Equation {

	private:
		KmerNb               _count;
		vector<unsigned int> _pathes;

	public:
		Equation () {}

		Equation (const KmerNb c): _count(c) {}

		KmerNb getNb() const {
			return _count;
		}

		const vector <unsigned int> &getPathes() const {
			return _pathes;
		}

		void addPath(unsigned int pathId) {
			_pathes.push_back(pathId);
		}

		KmerNb getScore(vector<KmerNb> &values) const {
			KmerNb score = 0;
			for (unsigned int i = 0; i < _pathes.size(); i++) {
				score += values[_pathes[i]];
			}
			return static_cast<KmerNb>(abs(_count - score));
		}
};


typedef map<int, Equation> Equations;

class EquationSystem {

	private:
		Equations            _equations;
		int                  _nbPathes;
		vector <KmerNb>      _values;
		vector<unsigned int> _sizes;
		KmerNb               _bestScore;

	public:
		EquationSystem (const int nbPathes);
		void addNode(const int nodeId, const KmerNb count);
		void addNodePath(const int pathId, const int nodeId);
		KmerNb getScore(vector <KmerNb> &values);
		KmerNb operator()(vector <KmerNb> &values);
		void fillSystem();
		void solveSimplex();
		
		KmerNb getValue(int index) const;

		friend ostream& operator<<(ostream& output, const EquationSystem& es);
};

#endif
