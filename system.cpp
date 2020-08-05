
#include "system.h"

System::System(double D, double k, double b,  int N,  double dx,
        std::vector<double> Pinit, std::vector<double> X)
: P(Pinit), J(N+1,0.), X(X), d(D), b(b), k(k), dx(dx), N(N), t(0)
{ }

void System::next_time(double dt)
{
    for(int i=1; i<N; ++i) {
        J[i] =  D((X[i-1]+X[i])/2.)*(P[i-1] - P[i])/dx;    
        J[i]+= -1*(A(X[i-1])*P[i-1] + A(X[i])*P[i])/2.;
    }
    // implement boundary
    J[0] = 0;
    J[N] = 0;
    //J[0] = -D*P[0]/dx;
    //J[N] = D*P[N-1]/dx;

    for(int i=0;i<N; ++i) {
        P[i] += dt*(J[i] - J[i+1])/dx;
    }
}





