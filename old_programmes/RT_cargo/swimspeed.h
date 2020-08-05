#ifndef GUARD_SWIMSPEED_H
#define GUARD_SWIMSPEED_H

class Swimspeed{
public:
    Swimspeed(double v0, double vp, double x0, double alpha)
            : v0(v0), vp(vp), x0(x0), alpha(alpha) {}


    //double v(double x) const { return v0 + vp*(x-x0); }
    double v(double x) const {
            if(x > x0 ) return v0 - vp*(x-x0);
            else return v0 + vp*(x-x0);
    }
    double get_alpha() const { return alpha; }

private:
    double v0, vp, x0, alpha;

};





#endif
