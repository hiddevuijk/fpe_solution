#ifndef GAURD_SYSTEM2D_H
#define GAURD_SYSTEM_2D_H

#include "xy.h"

#include <vector>
#include <math.h>

class System2d {
public:
    System2d(double D, double k, double b, int Nx, int Ny, double dx, double dy,
        std::vector<std::vector<double> > Pinit,
        std::vector<std::vector<double> > X);

    void next_time(double dt);

    std::vector<std::vector<double> > P;
    std::vector<std::vector<XY> > J;
    std::vector<std::vector<double> > X;

    double D(double x, double y)
        {return  d*(1+b*(x-2)*(x-2));}
    double A(double x, double y)
        { return k*x; }
private:
    double d;
    double b;
    double k;
    
    double dx;
    double dy;
    int Nx;
    int Ny;
    double t;
    


};

#endif
