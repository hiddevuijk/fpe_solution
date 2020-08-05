//////////////////////////////
//
// 2 dimensional vector class
// The +-*/ are defined as elementwise
// operations.
// 
//////////////////////////////

#ifndef GUARD_XY_H
#define GUARD_XY_H

#include <math.h>
#include <iostream>

class XY {
public:
	XY(double xx, double yy): x(xx), y(yy) {}
	XY(): x(0.), y(0.) {}

	// only two variable
	double x,y;

	// member functions	
	double length() const { return sqrt(x*x + y*y);}
	double length_sq() const {return x*x+y*y;}
	void normalize(double d=1.);
	void pbc(double L);
	XY operator+=(const XY& r);
	XY operator+=(const double& add);
	XY operator-=(const XY& r);
	XY operator-=(const double &sub);
	XY operator*=(const XY& r);
	XY operator*=(const double &mult);
	XY operator/=(const XY& r);
	XY operator/=(const double &div);
	
};

/////////////////////////
//
// NON MEMBER FUNCTIONS 
//
/////////////////////////


// out stream
std::ostream& operator<<(std::ostream &out, XY const& r);
// overload addition
XY operator+(const XY& l, const XY& r);
XY operator+(const XY& l, const double& add);
XY operator+(const double& add, const XY& r);
// overload subtraction
XY operator-(const XY& l, const XY& r);
XY operator-(const XY& l, const double& sub);
XY operator-(const double& sub, const XY& r);
// overload multiplication
XY operator*(const XY& l, const XY& r);
XY operator*(const XY& l, const double& mult);
XY operator*(const double& mult, const XY& r);
// overload division
XY operator/(const XY& l, const XY& r);
XY operator/(const XY& l, const double& div);


/////////////////////////
//
// NON MEMBER FUNCTIONS 
// IN NAME SPACE xy
//
/////////////////////////

namespace xy {
// distance between c1 and c2
double dist2(const XY &c1, const XY &c2);

// same as dist2, but with periodic boundary conditions
double dist2_pbc(const XY &c1, const XY &c2, double L);


// returns the distance squared
double dist2_sq(const XY &c1, const XY &c2);

// same as dist2_sq, but with periodic boundary conditions
double dist2_sq_pbc(const XY &c1, const XY &c2, double L);


// dot produc between c1 and c2
double dot(const XY &c1, const	XY &c2);

// the cross product, returns a double
double cross(const XY &c1, const XY &c2);


};

#endif
