
#include "system.h"

#include <iostream>
#include <vector>
#include <fstream>

using namespace std;

int main()
{


    int T = 100;
    int N = 50;


    double k = .01;
    double dt = 0.005;
    double dx = 1./N;
    double lambda = k*dt/(dx*dx);
    if(lambda >= 1./4) cerr << "fuck" << endl;


    vector<double> Pinit(N,0.);
    Pinit[N/2] = 0.5;
    Pinit[N/2-1] = 0.5;

    System system(lambda, N, T, Pinit);

    for(int ti=0; ti<T-1; ++ti)
        system.next_time(dt);

    for(int i=0;i<N; ++i) {
        for(int ti=0;ti<T; ti++) {
            cout << system.P[ti][i];
            if(ti < T-1) cout << '\t';
        }
        if( i < N-1) cout << '\n';
    }

        

    return 0;
}


