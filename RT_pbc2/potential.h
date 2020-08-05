#ifndef GUARD_POTENTIAL_H
#define GUARD_POTENTIAL_H

#include <math.h>

class Potential{
public:
    Potential(double k)
        : k(k){}
	
	//double U(double r) const { return k*r*r/2 ;}
	double F(double r) const { return -k*r; }
    
private:
    double k;
};

#endif
