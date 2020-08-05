#ifndef GUARD_SWIMSPEED_H
#define GUARD_SWIMSPEED_H

#include <math.h>

class Swimspeed{
public:
    Swimspeed(double v0, double vp, double x0, double alpha, double L)
            : v0(v0), vp(vp), x0(x0), alpha(alpha), L(L), L2(L/2)
			{}


    double v(double x, double y) const {
			double X = x+y;

			if( X > L2 ) X -= L;
			else if( X < -L2 ) X += L;

			if( X > L2 ) X -= L;
			else if( X < -L2 ) X += L;


			return v0 - vp*fabs(X-x0);
    }

    double get_alpha() const { return alpha; }

private:

    double v0, vp, x0, alpha;
	double L,L2;

};





#endif
