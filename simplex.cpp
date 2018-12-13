// Time-stamp: <simplex.cc   2 nov 2010 16:47:04> 

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <cstdlib>
#include <iostream>
#include "simplex.hpp"

// Constructs a dictionary
// Allocs arrays : coeffarray, baseofline and lineofbase
Simplex::Simplex(int nbConstraints, int nbVariables):
	_lin(nbConstraints + 1),               // Dimensions
	_col(nbVariables + nbConstraints + 1), 
	_coeffarray(_lin * _col, 0),           // Allocation and initialisation of the array of coefficients
	_baseofline(_lin, 0),                  // Allocation and initialisation of baseofline, no base variable at this time
	_lineofbase(_col, 0) { }               // Allocation and initialisation of lineofbase, no base variable

// Read accessor to a coefficient
const Rational& Simplex::getCoeff(unsigned int l, unsigned int c) const {
	if (l < _lin && c < _col) {
		return _coeffarray[l * _col + c];
	}
	else {
		cerr << "Error in Simplex::getCoeff: cell (" <<  l << ", " << c << ") is out of bound (" << _lin << ", " << _col << ")" << endl;
		return _coeffarray[0];
	}
}

// Write accessor to a coefficient
void Simplex::setCoeff(unsigned int l, unsigned int c, const Rational &r) {
	if (l < _lin && c < _col) {
		_coeffarray[l * _col + c] = r;
	}
	else {
		cerr << "Error in Simplex::setCoeff: cell (" <<  l << ", " << c << ") is out of bound (" << _lin << ", " << _col << ")" << endl;
	}
}

// Returns the base variable of a line
// Returns 
// -  0 if the line has no base variable
// -  j if the base variable of line i is x_j
int Simplex::baseOfLine(int l) {
	return _baseofline[l];
}

// Returns the line of a base variable
// Returns 
// -  0 if the variable x_j is not base variable
// -  i if x_j is the base variable of the line i 
int Simplex::lineOfBase(int c) {  
	return _lineofbase[c];
}

// Set the base variable of a base
void Simplex::setBase(int l, int c) {
	// Variable _baseofline[l] leaves the base
	if (_baseofline[l] != 0) {
		_lineofbase[_baseofline[l]] = 0;
	}
	// Variable c enters the base in line l
	_lineofbase[c] = l;
	_baseofline[l] = c;
}
    
  
// Return 
//  -  -1 if no variable can enter the basis
//  - the index of the first variable which can enter the basis
int Simplex::enter() {
	for (unsigned int c = 1; c < _col; c++) {
		if (getCoeff(0, c) > 0) {
			return c;
		}
	}
	return -1;			// No variable can enter the basis
}

// Return 
//  -  -1 if no variable can enter the basis
//  - the index of the first variable which the maximal positive coefficient
int Simplex::enterCoeff() {
	// Searching for the index 
	unsigned int max = 1;
	for (unsigned int c = 2; c < _col; c++) {
		if (getCoeff(0, c) > getCoeff(0, max)) {
			max = c;
		}
	}
	if (getCoeff(0, max) <= 0) {
		return -1;
	}
	return max;
}
	
  
// Computes the constraint on a variable given by a line
// The contraint is +infinity if the line number is -1
Rational Simplex::constraint(int l, int c) {
	if ((baseOfLine(l) == 0) || (getCoeff(l, c) <= 0)) {
		return Rational::infinity;
	}
	else {
		return getCoeff(l, 0) / getCoeff(l, c);
	}
}

// Constraint on a variable by a given line
// The contraint is +infinity if the line number is -1
Rational Simplex::constraint(int c) {
	Rational min = Rational::infinity;
	for (unsigned int l = 1; l < _lin; l++) {
		Rational val = constraint(l, c);
		if (val < min) {
			min = val;
		}
	}
	return min;
}

// Computes the line Most stringent constraint on a variable
// Return 
// -  -1 if there is no constraint
// -  the index of the lines which gives the most stringent constraint
int  Simplex::constline(int c) {
	Rational min = Rational::infinity;
	int minline = -1;		 // Line of the constraint
	for (unsigned int l = 1; l < _lin; l++) {
		Rational val = constraint(l, c);
		if (val < min) {
			min = val;
			minline = l;
		}
	}
	// No constraint
	if (minline == -1) {
		return -1;
	}
	return minline;
}
 
// Return 
//  -  -1 if no variable can enter the basis
//  - the index of the first variable which makes the increasing maximal
int Simplex::enterIncr() {
	// Searching for the index 
	Rational max = Rational::infinity.neg(); // -infinity
	int maxcol = -1;
	for (unsigned int c = 1; c < _col; c++) {
		if (getCoeff(0, c) > 0) {
			Rational val = getCoeff(0, c) * constraint(c);
			if (val > max) {
				max = val;
				maxcol = c;
			}
		}
	}
	if (maxcol == -1) {
		return -1;
	}
	return maxcol;
}

// Pivot operation     
void Simplex::pivot(unsigned int pl, unsigned int pc) {
	// Divides all element of the pivot raw by the pivot
	Rational pivot = getCoeff(pl, pc);
	for (unsigned int c = 0; c < _col; c++) {
		setCoeff(pl, c, getCoeff(pl, c) / pivot);
	}
	// Subtract the pivot row to all other raws
	for (unsigned int l = 0; l < _lin; l++) {
		if (l != pl) {
			Rational alpha = getCoeff(l, pc);
			for  (unsigned int c = 0; c < _col; c++) {
				setCoeff(l, c, getCoeff(l, c) - alpha * getCoeff(pl, c));
			}
		}
	}
	// Update base variables
	setBase(pl, pc);
}
	
Rational Simplex::solve() {
	while(true) {
		// Choice of the variable entering the basis
		int pc = enter();
		if (pc == -1) {
			return getCoeff(0, 0).neg();	// -a[0][0]
		}
		enterCoeff();
		enterIncr();
		// Choice of the variable leaving the basis
		int pl = constline(pc);
		if (pl == -1) {
			cerr << "Unbounded solution" << endl;
			return Rational::infinity; // +infinity
		}
		pivot(pl, pc);
	} 
}
 
// Returns true if the system has a solution
bool Simplex::init() {
	// Searching the minimal value of the b_i
	Rational min = Rational::infinity;
	int minline = -1;
	for (unsigned int l = 1; l < _lin; l++) {
		if (getCoeff(l, 0) < min) {
			min = getCoeff(l, 0);
			minline = l;
		}
	}
	// If all b_i are positive then the solution x_i = 0 satifies the equations
	if (min >= 0) {
		return true;
	}
	// Otherwise one solves an auxiliary problem
	Simplex aux = Simplex(_lin, _col-_lin);
	// Transfer this into aux
	// First line of aux : z = -x_0
	for (unsigned int c = 0; c < _col; c++) {
		aux.setCoeff(0, c, 0);
	}
	aux.setCoeff(0, _col, -1);
	// Last line of aux is the z-equation of this
	for (unsigned int c = 0; c < _col; c++) {
		aux.setCoeff(_lin, c, Rational(getCoeff(0, c)));
	}
	aux.setCoeff(_lin, _col, 0);
	// Other lines are as-is with -x_0
	for (unsigned int l = 1; l < _lin; l++) {
		for (unsigned int c = 0; c < _col; c++) {
			aux.setCoeff(l, c, Rational(getCoeff(l, c)));
		}
		aux.setCoeff(l, _col, -1);
	}
	// Setting bases
	for (unsigned int l = 1; l < _lin; l++) {
		aux.setBase(l, l + _col - _lin);
	}
	// Makes one pivot
	aux.pivot(minline, _col);
	// Solve the auxiliary problem
	if (aux.solve() < 0) {
		return false;	// Problem not feasible
	}
	// The auxiliary variable is the base
	if (aux.lineOfBase(_col) > 0) {
		// Searching a variable out of base to exchange them
		unsigned int c = 1;
		while (aux.lineOfBase(c++) != 0) // One variable out of base is found
			;
		aux.pivot(aux.lineOfBase(_col), c);
	}
	// Last line of aux is the z-equation of this
	for (unsigned int c = 0; c < _col; c++) {
		setCoeff(0, c, Rational(aux.getCoeff(_lin, c)));
	}
	// Other lines are copied back 
	for (unsigned int l = 1; l < _lin; l++) {
		for (unsigned int c = 0; c < _col; c++) {
			setCoeff(l, c, Rational(aux.getCoeff(l, c)));
		}
	}
	// Set base variables
	for (unsigned int l = 1; l < _lin; l++) {
		setBase(l, aux.baseOfLine(l));
	}
	return true;			// Feasible Simplex
}

Rational Simplex::simplex() {
	if (init()) {
		return solve();
	}
	else {
		cerr << "Unfeasable problem" << endl;
		return Rational::NaN;	// Problem not feasible
	}
}

void Simplex::addObjective(vector <int> &line) {
	if (line.size() != _col - _lin) {
		cerr << "Problem in simplex objective: expecting an array of size " << (_col - _lin) << ", got an array of size " << line.size() << endl;
	}
	for (unsigned int c = 1; c < _col - _lin + 1; c++) {
		setCoeff(0, c, Rational(line[c-1]));
	}
}

void Simplex::addLine(int lineNb, vector <int> &line, int b) {
	lineNb++;
	if (line.size() != _col - _lin) {
		cerr << "Problem in simplex constraint: expecting an array of size " << (_col - _lin) << ", got an array of size " << line.size() << endl;
	}
	for (unsigned int c = 1; c < _col - _lin + 1; c++) {
		setCoeff(lineNb, c, Rational(line[c-1]));
	}
	setCoeff(lineNb, 0, Rational(b));
	setCoeff(lineNb, lineNb + _col - _lin, 1);
}

void Simplex::setBase() {
	for (unsigned int l = 1; l < _lin; l++) {
		setBase(l, l + _col - _lin);
	}
}

float Simplex::getValue(int i) {
	i++;
	if (lineOfBase(i) == 0) {
		return 0;
	}
	return getCoeff(lineOfBase(i), 0);
}

// Scan Dictionary
istream& operator>>(istream& s, Simplex& simplex) {
	// Input the objective
	for (unsigned int c = 1; c < simplex._col - simplex._lin + 1; c++) {
		cout << "c[" << c << "] : ";
		s >> simplex._coeffarray[c];
	}
	// Scan equations
	for (unsigned int l = 1; l < simplex._lin; l++) {
		// Scan coefficients
		for (unsigned int c = 1; c < simplex._col - simplex._lin + 1; c++) {
			cout << "a[" << l << "," << c << "] : " ;
			s >> simplex._coeffarray[l * simplex._col + c];
		}
		cout << "b[" << l << "] : " ;
		s >> simplex._coeffarray[l * simplex._col];
		simplex.setCoeff(l, l+ simplex._col - simplex._lin, 1);
	}

	// Bases
	for (unsigned int l = 1; l < simplex._lin; l++) {
		simplex.setBase(l, l + simplex._col - simplex._lin);
	}
	return s;
}

// Print Dictionary
ostream& operator<<(ostream& s, Simplex& simplex) {
	// Print equations
	for (unsigned int l = 1; l < simplex._lin; l++) {
		if (simplex.baseOfLine(l) > 0)
			s << "x_" << simplex.baseOfLine(l) << " = ";
		else
			s << "      ";
		s << simplex.getCoeff(l, 0);
		for (unsigned int c = 1; c < simplex._col; c++) {
			if (simplex.lineOfBase(c) == 0) {
				if (simplex.getCoeff(l, c) < 0) {
					s << " + " << simplex.getCoeff(l, c).neg() << "x_" << c;
				}
				else if (simplex.getCoeff(l, c) == 0) {
					s << "       ";
				}
				else {
					s << " - " << simplex.getCoeff(l, c) << "x_" << c;
				}
			}
		}
		s << endl;
	}
	// Print objective
	s << "z   = " << simplex.getCoeff(0, 0).neg();
	for (unsigned int c = 1; c < simplex._col; c++) {
		if (simplex.lineOfBase(c) == 0) {
			if (simplex.getCoeff(0, c) > 0) {
				s << " + " << simplex.getCoeff(0, c)<< "x_" << c;
			}
			else if (simplex.getCoeff(0, c) == 0) {
				s << "       ";
			}
			else {
				s << " - " << simplex.getCoeff(0, c).neg()  << "x_" << c;
			}
		}
	}
	s << endl;
	return s;
}

