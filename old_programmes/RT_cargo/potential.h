#ifndef GUARD_POTENTIAL_H
#define GUARD_POTENTIAL_H

#include <math.h>

class Potential{
public:
    Potential(double k, double L = 0)
        : k(k), L2(L/2){}

    double U(double x, double y) const { return 0.5*k*(x - y)*(x - y); }
    double F(double x, double y) const { return k*(y - x); }
    //double F_pbc(double x, double y) const { return k*(y - x); }
    double F_pbc(double x, double y) const
	{ 
		return F(x,y);
		double d = fabs(x-y);
		if( d < L2) return F(x,y);
		else if( x > y ) return -1*F(x-L2,y);
		else return -1*F(x,y-L2);
		//if(d < L2) return k*(y-x);
		//else if( x > y ) {
		//	return -k*(x-y-L2);
		//} else {
		//	return k*(d-L2);
		//}
	}

private:
    double k;
	double L2;
};

#endif
