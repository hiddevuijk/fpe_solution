#ifndef GUARD_POTENTIAL_H
#define GUARD_POTENTIAL_H

class Potential{
public:
    Potential(double k, double kx0)
        : k(k), kx0(kx0) {}

    double U(double x) const { return 0.5*k*(x - kx0)*(x - kx0); }
    double F(double x) const { return k*(kx0 - x); }

private:
    double k, kx0;
};

#endif
