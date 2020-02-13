#include "system.h"

using namespace std;

System::System(double k, double v0, double vp, double x0, double d, double alpha, double D,
                double L, double dt, double dx, int N,
                vector<double> rhoInit, vector<double> sigmaInit,
                vector<double> PInit)
: k(k), v0(v0), vp(vp), x0(x0), d(d), alpha(alpha), D(D), L(L), dt(dt), dx(dx), N(N),
    x(N), rho(N), j0(N+1), sigma(N), j1(N+1), P(N), J(N+1)
{

    for(int i=0; i<N; ++i) {
        x[i] = 0.5*dx - 0.5*L + dx*i;
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

    // Bulk
    double v, xmid;
    //double vl, vr, xl, xr;
    for(int i=1; i<N; ++i) {
        xmid = (x[i] + x[i-1])/2.;
        //xl = x[i-1];
        //xr = x[i];
        v = v0 + vp*(xmid - x0);
        //vl = v0+ vp*(xl - x0);
        //vr = v0+ vp*(xr - x0);
        j0[i] = -d*( rho[i] - rho[i-1] )/dx;
        j0[i] += v*(sigma[i] + sigma[i-1])/2;
        //j0[i] += (vr*sigma[i] + vl*sigma[i-1])/2;

        j1[i] = -d*(sigma[i] - sigma[i-1])/dx;
        j1[i] += v*(rho[i] + rho[i-1])/2;
        //j1[i] += (vr*rho[i] + vl*rho[i-1])/2;

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

