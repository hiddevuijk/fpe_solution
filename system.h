#ifndef GAURD_SYSTEM_H
#define GAURD_SYSTEM_H

#include <vector>

class System {
public:
    System(double D, int N, double dx,  std::vector<double> Pinit);

    void next_time(double dt);

    std::vector<double> P;
    std::vector<double> J;
private:
    double D;
    double dx;
    int N;
    double t;
};

System::System(double D, int N,  double dx, std::vector<double> Pinit)
: D(D), N(N), dx(dx), P(N,0.), J(N+1,0.), t(0)
{
    for(int i=0; i<N; ++i)
        P[i] = Pinit[i];
}

void System::next_time(double dt)
{
    for(int i=1; i<N; ++i) {
        J[i] = D*(P[i-1] - P[i])/dx;    
    }
    // implement boundary
    // now open boundaries
    J[0] = 0;//D*P[0];
    J[N] = 0;//-D*P[N-1];

    for(int i=0;i<N; ++i) {
        P[i] += dt*(J[i] - J[i+1])/dx;
    }
}



#endif
