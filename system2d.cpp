
#include "system2d.h"

System2d::System2d(double D, double k, double b,  int Nx, int Ny,  double dx, double dy,
        std::vector<std::vector<double> > Pinit,
        std::vector<std::vector<double> > X)
: P(Pinit), J(Nx+1,std::vector<XY>(Ny+1) ), X(X), d(D), b(b), k(k),
    dx(dx), dy(dy), Nx(Nx), Ny(Ny), t(0)
{ }

void System2d::next_time(double dt)
{



}




