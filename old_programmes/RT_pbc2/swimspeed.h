#ifndef GUARD_SWIMSPEED_H
#define GUARD_SWIMSPEED_H

#include <math.h>

class Swimspeed{
public:
    Swimspeed(double v0, double vp, double x0, double alpha, double L)
            :pi2L(8*atan(1)/L), v0(v0), vp(vp), x0(x0), alpha(alpha), L(L), L2(L/2)
			{ }

    double v(double x) const {
       return v0+vp*sin(x*pi2L); 

    }
    //double v(double x) const {
	//		if( x > L2 ) x -= L;
	//		else if( x < -L2 ) x += L;

	//		if( x > L2 ) x -= L;
	//		else if( x < -L2 ) x += L;

    //        double pi2L = atan(1)*4/L;
	//		return v0 - vp*fabs(x-x0);
    //        //return fabs(v0*sin( (x-x0)*pi2L));
    //}

    double get_alpha() const { return alpha; }

private:
    double pi2L;
    double v0, vp, x0, alpha;
	double L,L2;

};





#endif
