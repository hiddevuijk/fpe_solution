#ifndef GUARD_SYSTEM_H
#define GUARD_SYSTEM_H

#include "swimspeed.h"
#include "potential.h"

#include <vector>

class System {
public:
    System( Swimspeed swimspeed, Potential potential,
            double d, double L, double dt, double dx, int N,
            std::vector<double> rhoInit,
            std::vector<double> sigmaInit);

    void next_time();

    const std::vector<double>& get_rho() const { return rho;}
    const std::vector<double>& get_sigma() const { return sigma;}
    const std::vector<double>& get_j0() const { return j0;}
    const std::vector<double>& get_j1() const { return j1;}

private:

    // PRIVATE MEMBER FUNCTIONS
    void next_flux();
    void next_prob();


    // SYSTEM PARAMETERS
    Swimspeed swimspeed;
    Potential potential;

    double d;       // diff. const. of RT
    double alpha;   // tumble rate


    double L;       // system size
    double dt;      // time step
    double dx;      // bin size
    int N;          // number of bins


    // THE STATE OF THE SYSTEM

    // bin centers
    std::vector<double> x;
    // prob. density of RT
    std::vector<double> rho;
    // zeroth moment of the flux
    std::vector<double> j0;
    // first moment of the prob. density
    std::vector<double> sigma;
    // first moment of the flux
    std::vector<double> j1;

};




#endif
