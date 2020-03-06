
#include "system.h"
#include "swimspeed.h"
#include "potential.h"

#include <iostream>
#include <fstream>
#include <vector>


using namespace std;

typedef vector<vector<double> > matrix;

int main()
{

    // potential
    double k = 0.;

    // swim/ tumble
    double v0 = 0.;
    double vp = 0.;
    double vx0 = 0;
    double alpha = 4.;

    // diffusion
    double gamma = 1.;
    double Gamma = 1.;
    double temp = 1.;

    

    double Lx = 10.;
    int Nx = 51;
    double dx = Lx/Nx;

    double Ly = 10.;
    int Ny = 51;
    double dy = Ly/Ny;

    double dt   = 0.00000001;
    double time = .05;
    int T = time/dt;
  

    vector<vector<double> > rInit(Nx, vector<double>(Ny,0));
    rInit[25][25] = 1./(Lx*Ly);
    vector<vector<double> > sInit(Nx, vector<double>(Ny,0));


    Swimspeed swimspeed(v0, vp, vx0, alpha);
    Potential potential(k);
    System system(swimspeed, potential, gamma, Gamma, temp, Nx, dx, Ny, dy, dt, rInit, sInit);

    double t = 0;
    int ti = 0;
    while( ti< T ) {
        system.next_time();
        t += dt;
        ti += 1;

        if( (ti % 100000) == 0) cout << T << '\t' << ti << endl;
    }

    ofstream out_r("r.dat");
    system.save_r(out_r);

    return 0;
}

