#ifndef GUARD_SYSTEM_H
#define GUARD_SYSTEM_H

#include <vector>

class System {
public:
    System(double k, double v0, double vp, double x0, double d, double alpha,
            double D, double L, double dt, double dx, int N,
            std::vector<double> rhoInit,
            std::vector<double> sigmaInit,
            std::vector<double> PInit);

    void next_time();

    const std::vector<double>& get_rho() const { return rho;}
    const std::vector<double>& get_sigma() const { return sigma;}
    const std::vector<double>& get_j0() const { return j0;}
    const std::vector<double>& get_j1() const { return j1;}

private:

    void next_flux();
    void next_prob();

    double k;       // spring const.
    double v0;      // swim speed
    double vp;      // slope swim speed
    double x0;      // zero of swim speed
    double d;       // diff. const. of RT
    double alpha;   // tumble rate

    double D;       // diff. const. of heavy particle

    double L;       // system size
    double dt;      // time step
    double dx;      // bin size
    int N;          // number of bins


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
    // prob. density of the heavy particle
    std::vector<double> P;
    // flux of the heavy particle
    std::vector<double> J;

};




#endif
