
#include "system.h"
#include "system2d.h"

#include <iostream>
#include <vector>
#include <fstream>

using namespace std;

int main()
{


    int T = 400000;
    int N = 51;

    double k = 0.05;
    double D = 0.1;
    double b = .05;
    double dt = 0.0001;

    double L = 10.;
    double dx = L/N;

    double lambda = D*dt/(dx*dx);
    if(lambda >= 1./2) cerr << " one" << endl;
    if(lambda >= 1./4) cerr << "FUCK: Lambda to large!" << endl;
    


    vector<double> Pinit(N,0.);
    Pinit[N/2 ] = 1.;

    vector<double> X(N,0.);
    for(int i=0;i<N;++i)
        X[i] = -0.5*L + dx*(i+0.5);

    System system(D, k,b, N, dx, Pinit, X);

    for(int ti=0; ti<T-1; ++ti)
        system.next_time(dt);

    for(int i=0;i<N; ++i) {
        cout <<  system.X[i] << '\t' << system.P[i];
        if( i < N-1) cout << '\n';
    }

        

    return 0;
}


