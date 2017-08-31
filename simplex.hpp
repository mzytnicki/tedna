// Time-stamp: <simplex.hh  29 Nov 2000 15:15:38>

#ifndef SIMPLEX_HPP
#define SIMPLEX_HPP

#include <vector>
#include "rational.hpp"

// Dictionary of a Linear problem
// The problem is represented using the tableau format
// cf. page 23 in "Linear Programming" of Vasek Chvatal 
//
// The slack variables are added to the original problem.
// If this original problem involves n variables and  m 
// equations,  the new problem has k = m+n variables and m 
// equations.
//
// The problem is stored in an array.  If the problem has
// m+n variables and m equations, the tableau has m+1 
// (from 0 to m) lines and m+n+1 columns (from 0 to m+n).
//
// Let A be the array where all coefficients are stored.
// The line 0 of A is devoted to the equation of z.
// The equation of z is written
// b_0 = c_1x_1 + .... + c_kx_k - z
// A[0][0] contains b_0
// A[0][j] contains c_j for j = 1..k
// The value of z is thus -A[0][0]
// The lines 1..m are devoted to equations.
// a_{i,1}x_1 + ... + a_{i,k}x_k = b_i for i = 1..m
// A[i][0] contains b_i for i = 1..m
// A[i][j] contains a_{i,j} for i = 1..m and j = 1..k
class Simplex {
public:
  Simplex(int nbConstraints, int nbVariables);		// Constructor
  Rational simplex();

  void addObjective(vector <int> &objective);
  void addLine(int lineNb, vector <int> &line, int b);
  void setBase();
  
  float getValue(int i);

  // Input and output
  friend istream& operator>>(istream& s, Simplex&); // Scanning
  friend ostream& operator<<(ostream& s, Simplex&); // Printing
private:
  // Handle entering and leaving variables
  int enter();
  int enterCoeff();
  Rational constraint(int l, int c);
  Rational constraint(int c);
  int constline(int c);
  int enterIncr();
  
  // Handle systems
  void pivot(unsigned int, unsigned int);
  Rational solve();
  bool init();
  
  // Handle Coefficients
  const Rational& getCoeff(unsigned int, unsigned int) const;	  // Read  accessor
  void  setCoeff(unsigned int, unsigned int, const Rational &); // Write accessor
  
  // Handle Base variables
  int baseOfLine(int);		// Return base variable of a line
  int lineOfBase(int);		// Return line of a base variable
  void setBase(int, int);	// Set the base variable of a line

  // Values
  unsigned int _lin;			// Number of lines
  unsigned int _col;			// Number of columns
  //int _lin;			// Number of lines
  //int _col;			// Number of columns

  // Array of coefficients
  // This one-dimentional array of size m x k contains the coefficients
  // The accessor coeff is used to access the coefficient (i,j)
  vector <Rational> _coeffarray;	
  
  // Array of base of each line
  // The value of baseofline[i] for i = 1..m is
  // -  0 if the line has no base variable
  // -  j if the base variable of line i is x_j
  // baseoline[0] is not used
  vector <int> _baseofline;		
  // Array of line of each base variable
  // The value of lineofbase[j] for j = 1..k is
  // -  0 if the variable x_j is not base variable
  // -  i if x_j is the base variable of the line i 
  // lineofbase[0] is not used
  vector <int> _lineofbase;
};

#endif 
