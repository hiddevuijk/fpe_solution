#include "xy.h"


/////////////////////////
//
// MEMBER FUNCTIONS
//
/////////////////////////

void XY::normalize(double d)
{	 double len = this->length();
	x *= d/len; y *= d/len;
}

void XY::pbc(double L)
	{ x -= L*floor(x/L); y -= L*floor(y/L); }

XY XY::operator+=(const XY& r)
	{ x += r.x; y += r.y; return *this; }

XY XY::operator+=(const double& add)
	{ x += add; y += add; return *this; }

XY XY::operator-=(const XY& r)
	{ x -= r.x; y -= r.y; return *this; }

XY XY::operator-=(const double& sub)
	{ x -= sub; y -= sub; return *this; }

XY XY::operator*=(const XY& r)
	{ x *= r.x; y *= r.y; return *this; }

XY XY::operator*=(const double& mult)
	{ x *= mult; y *= mult; return *this; }

XY XY::operator/=(const XY& r)
	{ x /= r.x; y /= r.y; return *this; }

XY XY::operator/=(const double& div)
	{ x /= div; y /= div; return *this; }


	
/////////////////////////
//
// NON MEMBER FUNCTIONS 
//
/////////////////////////


// out stream
std::ostream& operator<<(std::ostream &out, XY const& r)
{ out << r.x << '\t' << r.y ; return out; }

// overload addition
XY operator+(const XY& l, const XY& r)
{ return XY(l.x+r.x, l.y+r.y); }

XY operator+(const XY& l, const double& add)
{ return XY(l.x+add, l.y+add); }

XY operator+(const double& add, const XY& r)
{ return XY(r.x+add, r.y+add); }

// overload subtraction
XY operator-(const XY& l, const XY& r)
{ return XY(l.x-r.x, l.y-r.y); }

XY operator-(const XY& l, const double& sub)
{ return XY(l.x-sub, l.y-sub); }

XY operator-(const double& sub, const XY& r)
{ return XY(sub-r.x,sub-r.y); }

// overload multiplication
XY operator*(const XY& l, const XY& r)
{ return XY(l.x*r.x, l.y*r.y); }

XY operator*(const XY& l, const double& mult)
{ return XY(l.x*mult, l.y*mult); }

XY operator*(const double& mult, const XY& r)
{ return XY(mult*r.x,mult*r.y); }


// overload division
XY operator/(const XY& l, const XY& r)
{ return XY(l.x/r.x, l.y/r.y); }

XY operator/(const XY& l, const double& div)
{ return XY(l.x/div, l.y/div); }




/////////////////////////
//
// NON MEMBER FUNCTIONS 
// IN NAME SPACE xy
//
/////////////////////////

// distance between c1 and c2
double xy::dist2(const XY &c1, const XY &c2)
{ return sqrt( (c1.x - c2.x)*(c1.x - c2.x) +
			(c1.y - c2.y)*(c1.y - c2.y) ); }

// same as dist2, but with periodic boundary conditions
double xy::dist2_pbc(const XY &c1, const XY &c2, double L)
{
	XY d;
	d.x = c1.x - c2.x;
	d.y = c1.y - c2.y;
	d.x -= L*round(d.x/L);
	d.y -= L*round(d.y/L);
	return sqrt(d.x*d.x+d.y*d.y);
}


// returns the distance squared
double xy::dist2_sq(const XY &c1, const XY &c2)
{ return  (c1.x - c2.x)*(c1.x - c2.x) +
			(c1.y - c2.y)*(c1.y - c2.y);
}

// same as dist2_sq, but with periodic boundary conditions
double xy::dist2_sq_pbc(const XY &c1, const XY &c2, double L)
{
	XY d;
	d.x = c1.x - c2.x;
	d.y = c1.y - c2.y; 
	d.x -= L*round(d.x/L);
	d.y -= L*round(d.y/L); 
	return d.x*d.x + d.y*d.y;
}


// dot produc between c1 and c2
double xy::dot(const XY &c1, const	XY &c2) 
{ return c1.x*c2.x + c1.y*c2.y; }

// the cross product, returns a double
double xy::cross(const XY &c1, const XY &c2) 
{ return c1.x*c2.y - c1.y*c2.x; }


