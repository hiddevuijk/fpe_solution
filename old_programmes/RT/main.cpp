
#include "system.h"
#include "swimspeed.h"
#include "potential.h"

#include <iostream>
#include <vector>


using namespace std;

int main()
{


    double k = 1.;
    double kx0 = 0;
    double v0 = 5.;
    double vp = 0.;
    double vx0 = -2.5;
    double d = 1.;
    double alpha = 40.;


    double L = 10.;
    double dt = 0.00001;
    int N = 101;
    double dx = L/N;

    double time = 10;
    int T = time/dt;


    vector<double> rhoInit(N, 1./L);
    vector<double> sigmaInit(N, 0.0);

    Swimspeed swimspeed(v0, vp, vx0, alpha);
    Potential potential(k, kx0);
    System system(swimspeed, potential, d, L, dt, dx, N, rhoInit, sigmaInit);

    for(int i=0; i<T; ++i) 
        system.next_time();

    for(int i=0;i<N; ++i) {
        cout <<  0.5*dx -L/2. + i*dx << '\t'
             << system.get_rho()[i] << '\t'
             << system.get_sigma()[i] << '\t'
             << system.get_j0()[i] << '\t'
             << system.get_j1()[i] << '\n';
    }


    return 0;
}

