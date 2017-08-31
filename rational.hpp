// Time-stamp: <rational.hh   2 nov 2010 16:41:27>

#ifndef RATIONAL_HPP
#define RATIONAL_HPP

#include <iostream>
using namespace std;

// Class to implement rational integers
class Rational {
public:
  Rational(int = 0, int = 1);	             // Constructor from integers
  Rational(const Rational&);	             // Constructor by copy
  Rational& operator=(const Rational&);      // Affectation
  Rational& operator=(int n);                // Affectation
  Rational neg () const;		             // Opposite (unary -)
  bool operator==(const Rational&) const;    // Equality
  bool operator==(int) const;                // Equality
  bool operator<(const Rational&) const;     // Order (strict)
  bool operator<(int) const;	             // Order (strict)
  bool operator<=(const Rational&) const;    // Order
  bool operator<=(int) const;                // Order
  bool operator>(const Rational&) const;     // Order (strict)
  bool operator>(int) const;	             // Order (strict)
  bool operator>=(const Rational&) const;    // Order
  bool operator>=(int) const;                // Order
  Rational operator+(const Rational&) const; // Addition
  Rational operator-(const Rational&) const; // Difference
  Rational operator*(const Rational&) const; // Multiplication
  Rational operator/(const Rational&) const; // Division

  operator float() const;     // export to float
  
  static const Rational NaN;           // Constant for NaN
  static const Rational infinity;      // Constant for infinity

  friend istream& operator>>(istream& s, Rational&); // Scanning
  friend ostream& operator<<(ostream& s, const Rational&); // Printing
private:
  int num;			// Numerator
  int den;			// Denominator
  void canonical();		// Make the representation canonical
};

#endif 
