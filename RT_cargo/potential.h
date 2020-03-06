#ifndef GUARD_POTENTIAL_H
#define GUARD_POTENTIAL_H

class Potential{
public:
    Potential(double k)
        : k(k){}

    double U(double x, double y) const { return 0.5*k*(x - y)*(x - y); }
    double F(double x, double y) const { return k*(y - x); }

private:
    double k;
};

#endif
