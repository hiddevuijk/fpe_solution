
#include "system.h"
#include "swimspeed.h"
#include "potential.h"

#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>


using namespace std;

typedef vector<vector<double> > matrix;

int main()
{

    // system size
    double Lx = 10.;
    int Nx = 101;
    double dx = Lx/Nx;

    double Ly = 10.;
    int Ny = 101;
    double dy = Ly/Ny;

    // potential
    double k = 1.;

    // swim/ tumble
    double v0 = 10.;
    double vp = 2*v0/Lx;
    double vx0 = 0.;

    //double v0 = 0.;
    //double vp = 1;
    //double vx0 = -Lx/2.;

    double alpha = 10.;

    // diffusion
    double gamma = .5;
    double Gamma = 10.;
    double temp = 1.;

    // integration time 
    double dt   = 0.001;
    double time = 100.;
    int T = time/dt;

	cout << "Neumann criterion: \n";	
	cout << "x: \t" << temp*dt/(gamma*dx*dx) << "\n";	
	cout << "y: \t" << temp*dt/(Gamma*dy*dy) << "\n";	
  

    vector<vector<double> > rInit(Nx, vector<double>(Ny,0));
    rInit[int(Nx/2.)][int(Ny/2.)] = 1.;
    vector<vector<double> > sInit(Nx, vector<double>(Ny,0));


    Swimspeed swimspeed(v0, vp, vx0, alpha);
    Potential potential(k);
    System system(swimspeed, potential, gamma, Gamma, temp, Nx, dx, Ny, dy, dt, rInit, sInit);

    double t = 0;
    int ti = 0;
    int ni = 0;
    while( ti< T ) {

        if( (ti % int(T/10)) == 0) {
            
            cout << T << '\t' << ti << endl;

            stringstream sr;
            sr << "r" << ni << ".dat";
            ofstream out_r(sr.str());
            system.save_r(out_r);

            stringstream ss;
            ss << "s" << ni << ".dat";
            ofstream out_s(ss.str());
            system.save_s(out_s);

            ni += 1;

        }

        system.next_time();
        t += dt;
        ti += 1;
    }

    ofstream out_r("r.dat");
    system.save_r(out_r);

    ofstream out_s("s.dat");
    system.save_s(out_s);

    ofstream out_x("x.dat");
    system.save_x(out_x);

    ofstream out_y("y.dat");
    system.save_y(out_y);
    return 0;
}

