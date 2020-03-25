#ifndef GUARD_SYSTEM_H
#define GUARD_SYSTEM_H

#include "swimspeed.h"
#include "potential.h"

#include <vector>
#include <fstream>

class System {
public:
    System( Swimspeed swimspeed, Potential potential,
            double gamma, double Gamma, double temp, bool periodic,
			int Nx, double dx, int Ny, double dy, double dt,
            std::vector<std::vector<double> > rInit,
            std::vector<std::vector<double> > sInit);

	// increment one time step
    void next_time();

	// get state of the system
    const std::vector<std::vector<double> >& get_r() const { return r;}
    const std::vector<std::vector<double> >& get_s() const { return s;}
    const std::vector<std::vector<double> >& get_jx0() const { return jx0;}
    const std::vector<std::vector<double> >& get_jy0() const { return jy0;}
    const std::vector<std::vector<double> >& get_jx1() const { return jx1;}
    const std::vector<std::vector<double> >& get_jy1() const { return jy1;}

	// write the state to file
    void save_r(std::ofstream& out) const;
    void save_s(std::ofstream& out) const;
    void save_x(std::ofstream& out) const;
    void save_y(std::ofstream& out) const;

private:

    // PRIVATE MEMBER FUNCTIONS
    void next_flux();
    void next_flux_pbc();
    void next_prob();
    void next_prob_pbc();


    // SYSTEM PARAMETERS
    Swimspeed swimspeed;
    Potential potential;

    double gamma;
    double Gamma;
    double d, D;       // diff. const. of RT
    double temp;

    bool periodic;  // periodic BC
    double Lx, Ly;       // system size
    int Nx;          // number of bins
    double dx;
    int Ny;
    double dy;      // bin size
    double dt;      // time step


    // THE STATE OF THE SYSTEM

    // bin centers
    std::vector<double> x;
    std::vector<double> y;
    // prob. density of RT
    std::vector<std::vector<double> > r;
    // zeroth moment of the flux
    std::vector<std::vector<double> > jx0;
    std::vector<std::vector<double> > jy0;
    // first moment of the prob. density
    std::vector<std::vector<double> > s;
    // first moment of the flux
    std::vector<std::vector<double> > jx1;
    std::vector<std::vector<double> > jy1;


};




#endif
