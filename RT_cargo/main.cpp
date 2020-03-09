
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
    double k = 1.;

    // swim/ tumble
    double v0 = 0.;
    double vp = 1.;
    double vx0 = 0.;
    double alpha = 10.;

    // diffusion
    double gamma = 5.;
    double Gamma = .5;
    double temp = 1.;

    

    double Lx = 10.;
    int Nx = 50;
    double dx = Lx/Nx;

    double Ly = 10.;
    int Ny = 51;
    double dy = Ly/Ny;

    double dt   = 0.00001;
    double time = 10.;
    int T = time/dt;
  

    vector<vector<double> > rInit(Nx, vector<double>(Ny,0));
    rInit[25][25] = 1.;
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

        if( (ti % int(T/10)) == 0) cout << T << '\t' << ti << endl;
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

