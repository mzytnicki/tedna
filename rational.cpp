#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

// Time-stamp: <rational.cc   2 nov 2010 16:34:59> 

#include "rational.hpp"
#include <cstdlib>

// Constants
// Constant for NaN
const Rational Rational::NaN = Rational(0, 0);
// Constant for infinity
const Rational Rational::infinity = Rational(1, 0);

// Computes the sign of an integer
//  -  -1   if  n < 0
//  -   0   if  n = 0
//  -   1   if  n > 0
int sign(int n) {
	if (n < 0) 
		return -1;
	else if (n == 0)
		return 0;
	else
		return 1;
}
    
// Computes the GCD of the two integers.
// m and are assumed to be non-negative.
int gcd(int m, int n) {
	int r;			// Temporary variable

	// Exchange m and n if m < n
	if (m < n) {
		r = n;  n = m; m = r;
	}
	// It can be assumed that m >= n
	while (n > 0) {
		r = m % n;
		m = n;
		n = r;
	}
	return m;
}

// Makes the representation canonical
// The canonical representations are the following
//   -   1/0  +infinity
//   -  -1/0  -infinity
//   -   0/0  NaN
//   -   0/1  0
//   -   n/d  with n != 0, d > 0 and n and d coprime
void Rational::canonical() {
	if (den == 0) {
		// Special numbers : +/-infinity and NaN
		num = sign(num);
		return;
	}
	if (num == 0) {
		// 0
		den = 1;
		return;
	}
	int sign = 1;			// Sign of the rational
	// Makes num positive and updates sign
	if (num < 0) {
		num = -num;
		sign = -1;
	}
	// Makes den positive and updates sign
	if (den < 0) {
		den = -den;
		sign = -sign;
	}
	int d = gcd(num, den);	// Gcd of num and den
	num = sign*num/d;
	den = den/d;
}

// Constructs a Rational from the numerator and denominator.
// This constructor is also the default constructor since
// both n and d have default values.  This default constructor
// is needed to construt arrays of Rationals.
Rational::Rational(int n, int d) {
	num = n;
	den = d;
	canonical();
}

// Constructs a Rational from another Rational
// It is just the copy of the numerator and denominator.
Rational::Rational(const Rational& r) {
	num = r.num;
	den = r.den;
}

// Affectation 
// It is just the copy of the numerator and denominator.
Rational& Rational::operator=(const Rational& r) {
	num = r.num;
	den = r.den;
	return *this;
}

// Affectation 
// It is just the copy of the numerator and the denominator 
// is set to 1.
Rational& Rational::operator=(int n) {
	num = n;
	den = 1;
	return *this;
}

// Opposite
Rational Rational::neg() const {
	return Rational(-num, den);
}

// Equality
bool Rational::operator==(const Rational& r) const {
	return num == r.num && den == r.den;
}

// Equality with an integer
bool Rational::operator==(int n) const {
	return num == n && den == 1;
}

// Orders
// All orders are defined from the order <
// 	     NaN    -infty    +infty    p'/q'
// NaN        F        F         F        F
// -infty     F        F         T        T
// +infty     F        F         F        F
// p/q        F        F         T     pq'<p'q
bool Rational::operator<(const Rational& r) const {
	if (den == 0 && r.den == 0) {
		// Both are special
		if (num == 0 || r.num == 0)
			// One of them is NaN
			return 0;
		else
			// Both are -infty or +infty
			return num < r.num;
	} else
		// Other cases
		return num*r.den < den*r.num;
}

bool Rational::operator<(int n) const {
	return *this < Rational(n);
}

// Equality or Order 
bool Rational::operator<=(const Rational& r) const {
	return *this == r || *this < r;
}

bool Rational::operator<=(int n) const {
	return *this <= Rational(n);
}

// Reversed order
bool Rational::operator>(const Rational& r) const {
	return r < *this;
}

bool Rational::operator>(int n) const {
	return Rational(n) < *this;
}

// Equality or reversed order
bool Rational::operator>=(const Rational& r) const {
	return *this == r || r < *this;
}

bool Rational::operator>=(int n) const {
	return Rational(n) <= *this;
}

// Addition
// 	     NaN    -infty    +infty    p'/q'
// NaN       NaN      NaN       NaN      NaN
// -infty    NaN    -infty      NaN    -infty
// +infty    NaN      NaN     +infty   +infty
// p/q       NaN    -infty    -infty  pq'+p'q/dd'
Rational Rational::operator+(const Rational& r) const {
	if (den == 0 || r.den == 0) {
		if (den == 0 && r.den == 0) {
			// Both are special
			if (num == 0 || r.num == 0)
				// One is NaN
				return Rational(0, 0);	// NaN
			else
				// None is NaN
				return Rational(num + r.num, 0);
		} else 
			return Rational(num*r.den+den*r.num, 0);
	}
	else {
		int d = gcd(den, r.den);	// Gcd of denominators
		int dens = den/d;		// Simplified denominator of this
		int rdes = r.den/d;		// Simplified denominator of r
		return Rational(num*rdes + dens * r.num, den * rdes);
	}
}

// Difference
Rational Rational::operator-(const Rational& r) const {
	if (den == 0 || r.den == 0) {
		if (den == 0 && r.den == 0) {
			// Both are special
			if (num == 0 || r.num == 0)
				// One is NaN
				return Rational(0, 0);	// NaN
			else
				// None is NaN
				return Rational(num - r.num, 0);
		}
		else 
			return Rational(num*r.den - den*r.num, 0);
	} else {
		int d = gcd(den, r.den);	// Gcd of denominators
		int tdes = den/d;		// Simplified denominator of this
		int rdes = r.den/d;		// Simplified denominator of r
		return Rational(num*rdes - tdes * r.num, den * rdes);
	}
}

// Multiplication
// 	     NaN    -infty    +infty      0      p'/q'
// NaN       NaN      NaN       NaN      NaN      NaN
// -infty    NaN    -infty      NaN      NaN    -infty     
// +infty    NaN      NaN     +infty     NaN    +infty
//  0        NaN      NaN       NaN       0        0
// p/q       NaN    -infty    -infty      0    pq'+p'q/dd'
Rational Rational::operator*(const Rational& r) const {
	return Rational(num*r.num, den*r.den);
}

// Division
Rational Rational::operator/(const Rational& r) const {
	return Rational(num*r.den, den*r.num);
}

// Export to float
Rational::operator float() const {
	return (num/den);
}

// Output a Rational
istream& operator>>(istream& s, Rational& r) {
	char buffer[256];
	char* denstr;

	// Input a string
	s >> buffer;

	// Search for a '/' until end of string 
	for (denstr = buffer; *denstr != '\0' && *denstr != '/'; denstr++);

	// If the string contains '/', scan denominator
	if (*denstr == '/') {
		*denstr++ = '\0';		// Marks end of numerator
		r.den = atoi(denstr);	// Scan denominator
	} else { 
		r.den = 1;			// Otherwise denominator is 1
	}
	r.num = atoi(buffer);
	r.canonical();

	return s;
}

// Output a Rational
ostream& operator<<(ostream& s, const Rational& r) {
	switch (r.den) {
		case 0:
			switch (r.num) {
				case 1:
					s << "+infty";
					return s;
				case -1:      
					s << "-infty";
					return s;
				case 0:      
					s << "NaN";
					return s;
				default:
					s << "Non canonical Rational";
					return s;
			}
		case 1:
			s << r.num;
			return s;
		default:
			s << r.num << "/" << r.den;
			return s;
	}
}
