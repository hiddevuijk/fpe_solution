
#include "system.h"

#include <iostream>
#include <vector>


using namespace std;

int main()
{


    double k = 1.;
    double v = 1.;
    double d = 5.;
    double alpha = 1.;

    double D = 1.;

    double L = 5.;
    double dt = 0.001;
    double dx = 0.1;

    int N = L/dx;

    vector<double> rhoInit(N,0.0);
    vector<double> sigmaInit(N,0.0);
    vector<double> PInit(N,0);

    rhoInit[25] = 1./dx;

    System system(k,v,d,alpha,D,L,dt,dx,N, rhoInit, sigmaInit, PInit);

    for(int i=0; i<5000000; ++i) 
        system.next_time();

    for(int i=0;i<N;++i) {
        cout << system.rho[i] << '\t'
            << system.sigma[i] << endl;
    }



    return 0;
}

