#ifndef GAURD_SYSTEM_H
#define GAURD_SYSTEM_H

#include <vector>
#include <math.h>

class System {
public:
    System(double D, double k, double b, int N, double dx,  
        std::vector<double> Pinit, std::vector<double> X);

    void next_time(double dt);

    std::vector<double> P;
    std::vector<double> J;
    std::vector<double> X;

    double D(double x)
        {return  d*(1+b*(x-2)*(x-2));}
    double A(double x)
        { return k*x; }
private:
    double d;
    double b;
    double k;
    
    double dx;
    int N;
    double t;
    


};

#endif
