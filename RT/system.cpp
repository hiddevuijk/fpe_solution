#include "system.h"

using namespace std;

System::System(double k, double v, double d, double alpha, double D,
                double L, double dt, double dx, int N,
                vector<double> rhoInit, vector<double> sigmaInit,
                vector<double> PInit)
: k(k), v(v), d(d), alpha(alpha), D(D), L(L), dt(dt), dx(dx), N(N),
    rho(N), j0(N+1), sigma(N), j1(N+1), P(N), J(N+1)
{

    for(int i=0; i<N; ++i) {
        rho[i] = rhoInit[i];
        sigma[i] = sigmaInit[i];
        P[i] = PInit[i];
    }

}


void System::next_time()
{
    next_flux();
    next_prob();
}


void System::next_flux()
{
    // BC
    j0[0] = 0;     
    j0[N] = 0;

    j1[0] = 0;     
    j1[N] = 0;


    for(int i=1; i<N; ++i) {
        // Bulk
        j0[i] = -d*( rho[i] - rho[i-1] )/dx;
        j0[i] += v*(sigma[i] + sigma[i-1])/(2*dx);

        j1[i] = -d*(sigma[i] - sigma[i-1])/dx;
        j1[i] += v*(rho[i] + rho[i-1])/(2*dx);

    }

}

void System::next_prob()
{

    for(int i=0;i<N; ++i) {
        rho[i] += dt*j0[i]/dx;
        rho[i] -= dt*j0[i+1]/dx;

        sigma[i] -= dt*alpha*sigma[i];
        sigma[i] += dt*j1[i]/dx;
        sigma[i] -= dt*j1[i+1]/dx;
    }

}

