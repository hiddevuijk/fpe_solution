#ifndef GAURD_SYSTEM_H
#define GAURD_SYSTEM_H

#include <vector>

class System {
public:
    System(double lambda, int N, int T, std::vector<double> Pinit);

    void next_time(double dt);

    std::vector<std::vector<double> > P;
private:
    double lambda;
    int N;
    int T;
    int ti;
};

System::System(double lambda, int N, int T, std::vector<double> Pinit)
: lambda(lambda), N(N),T(T), P(T, std::vector<double>(N,0.) ), ti(0)
{
    for(int i=0; i<N; ++i)
        P[0][i] = Pinit[i];


}

void System::next_time(double dt)
{
    for(int i=1; i<N-1; ++i) {
        P[ti+1][i] = P[ti][i] +
                 lambda*( P[ti][i-1] - 2*P[ti][i] + P[ti][i+1] );
    }
    ++ti;

}



#endif
