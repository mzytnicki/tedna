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

#include <cmath>
#include "globals.hpp"
#include "equations.hpp"
#include "simplex.hpp"


EquationSystem::EquationSystem(const int nbPathes): _nbPathes(nbPathes), _sizes(_nbPathes, 0) {}

void EquationSystem::addNode(const int nodeId, const KmerNb count) {
	_equations[nodeId] = Equation(count);
}

void EquationSystem::addNodePath(const int pathId, const int nodeId) {
	_equations[nodeId].addPath(pathId);
	_sizes[pathId]++;
}

KmerNb EquationSystem::getScore(vector <KmerNb> &values) {
	KmerNb score = 0;
	for (Equations::const_iterator it = _equations.begin(); it != _equations.end(); ++it) {
		score += it->second.getScore(values);
	}
	return score;
}

KmerNb EquationSystem::operator()(vector <KmerNb> &values) {
	return getScore(values);
}

void EquationSystem::solveSimplex() {
	int nbPathes    = _nbPathes;
	int nbNodes     = _equations.size();
	int nbVariables = nbNodes + nbPathes;
	int nbEquations = 2 * nbNodes;
	Simplex simplex = Simplex(nbEquations, nbVariables);
	/*
	vector <int> objective(nbPathes, 0);
	for (int i = 0; i < _nbPathes; i++) {
		objective[i] = _sizes[i] * _sizes[i];
	}
	*/
	vector <int> objective(nbPathes, 0);
	vector <int> objective2(nbNodes, -1);
	objective.insert(objective.end(), objective2.begin(), objective2.end());
	simplex.addObjective(objective);
	int i = 0;
	for (Equations::const_iterator it = _equations.begin(); it != _equations.end(); ++it, i++) {
		vector <int> line(nbPathes, 0);
		vector <int> line2(nbPathes, 0);
		vector <int> line3(nbNodes, 0);
		const vector <unsigned int> &pathes = it->second.getPathes();
		for (unsigned int j = 0; j < pathes.size(); j++) {
			line[pathes[j]] = 1;
			line2[pathes[j]] = -1;
		}
		line3[i] = -1;
		line.insert(line.end(), line3.begin(), line3.end());
		line2.insert(line2.end(), line3.begin(), line3.end());
		simplex.addLine(2 * i, line, it->second.getNb());
		simplex.addLine(2 * i + 1, line2, - it->second.getNb());
	}
	simplex.setBase();
	//cout << "Simplex in: " << endl << simplex << endl;
	simplex.simplex();
	//cout << "Simplex out: " << endl << simplex << endl;
	//for (int i = 0; i < nbPathes+nbNodes; i++) {
	//	cout << simplex.getValue(i) << " ";
	//}
	//cout << endl;
	for (int i = 0; i < nbPathes; i++) {
		_values.push_back(simplex.getValue(i));
	}
}


KmerNb EquationSystem::getValue(int index) const {
	return _values[index];
}

ostream& operator<<(ostream& output, const EquationSystem& es) {
	for (unsigned int i = 0; i < es._values.size(); i++) {
		output << es._values[i] << "  ";
	}
	output << endl;
	return output;
}
