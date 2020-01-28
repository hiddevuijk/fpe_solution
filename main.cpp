
#include "system.h"

#include <iostream>
#include <vector>
#include <fstream>

using namespace std;

int main()
{


    int T = 50000;
    int N = 50;


    double D = 0.1;
    double dt = 0.005;

    double L = 10.;
    double dx = L/N;

    double lambda = D*dt/(dx*dx);
    if(lambda >= 1./4) cerr << "FUCK: Lambda to large!" << endl;
    


    vector<double> Pinit(N,0.);
    Pinit[N/2 - 1] = 1.;

    System system(D, N, dx, Pinit);

    for(int ti=0; ti<T-1; ++ti)
        system.next_time(dt);

    for(int i=0;i<N; ++i) {
        cout <<  (i+0.5)*dx << '\t' << system.P[i];
        if( i < N-1) cout << '\n';
    }

        

    return 0;
}


