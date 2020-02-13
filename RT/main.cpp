
#include "system.h"

#include <iostream>
#include <vector>


using namespace std;

int main()
{


    double k = 1.;
    double v0 = 0.;
    double vp = 1.;
    double x0 = -2.5;
    double d = 1.;
    double alpha = .1;

    double D = 1.;

    double L = 5.;
    double dt = 0.00001;
    int N = 501;
    double dx = L/N;

    double time = 10;
    int T = time/dt;


    vector<double> rhoInit(N,1./L);
    vector<double> sigmaInit(N,0.0);
    vector<double> PInit(N,0);


    System system(k,v0, vp, x0, d,alpha,D,L,dt,dx,N, rhoInit, sigmaInit, PInit);

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

