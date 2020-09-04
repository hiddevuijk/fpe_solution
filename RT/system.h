#ifndef GUARD_SYSTEM_H
#define GUARD_SYSTEM_H

#include "swimspeed.h"
#include "potential.h"

#include <vector>
#include <fstream>
#include <string>

class System {
public:
    System( Swimspeed swimspeed, Potential potential,
            double gamma, double Gamma, double temp,
			int Nx, double dx, int Ny, double dy,
            const std::vector<std::vector<double> >& rInit,
            const std::vector<std::vector<double> >& sInit);

	// increment one time step
    void next_time(double dt);
    double next_time_check_error(double dt);

    void set_init();
    void read_init(std::string);

	// get state of the system
    std::vector<std::vector<double> > get_r() const { return r;}
	std::vector<double> get_rx() const;
	std::vector<double> get_ry() const;
    std::vector<std::vector<double> > get_s() const { return s;}
    std::vector<std::vector<double> > get_jx0() const { return jx0;}
    std::vector<std::vector<double> > get_jy0() const { return jy0;}
    std::vector<std::vector<double> > get_jx1() const { return jx1;}
    std::vector<std::vector<double> > get_jy1() const { return jy1;}

	double get_r(int xi, int yi) const { return r[xi][yi]; }
	double get_s(int xi, int yi) const { return s[xi][yi]; }



	// write the state to file
	void save_state(std::string path, std::string name) const;
	void save_state_1d(std::string path, std::string name) const;

	void save_flux(std::string path, std::string name) const;
	void save_xy(std::string path, std::string name) const;

    void save_r(std::ofstream& out) const;
    void save_s(std::ofstream& out) const;
    void save_x(std::ofstream& out) const;
    void save_y(std::ofstream& out) const;
	void save_jx0(std::ofstream& out) const;
	void save_jy0(std::ofstream& out) const;
	void save_jx1(std::ofstream& out) const;
	void save_jy1(std::ofstream& out) const;

	// parameters
	int get_Nx() const {return Nx;}
	int get_Ny() const {return Ny;}
	double get_dx() const {return dx;}
	double get_dy() const {return dy;}

private:

    // PRIVATE MEMBER FUNCTIONS
    void next_flux();
    void next_prob(double dt);


    // SYSTEM PARAMETERS
    Swimspeed swimspeed;
    Potential potential;

    double gamma;
    double Gamma;
    double q;
    double d, D;       // diff. const. of RT
    double temp;

    int Nx;          // number of bins
    double dx;
    int Ny;
    double dy;      // bin size


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

double steady_state_error( const std::vector<std::vector<double> >& r,
				   const std::vector<std::vector<double> >& s,
				   const System& system);


bool steady_state( const std::vector<std::vector<double> >& r,
				   const std::vector<std::vector<double> >& s,
				   const System& system, double epsilon);

#endif
