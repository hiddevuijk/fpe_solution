#include "system.h"

using namespace std;

System::System( Swimspeed swimspeed, Potential potential,
                double gamma, double Gamma, double temp, int Nx, double dx, int Ny, double dy, double dt,
                vector<vector<double> > rInit, vector<vector<double> > sInit)
: swimspeed(swimspeed), potential(potential),
    gamma(gamma), Gamma(Gamma), temp(temp),
    Nx(Nx), dx(dx), Ny(Ny), dy(dy), dt(dt),
    x(Nx), y(Ny),
    r( Nx, std::vector<double>(Ny) ),
    jx0( Nx+1, std::vector<double>(Ny+1) ),
    jy0( Nx+1, std::vector<double>(Ny+1) ),
    s( Nx, std::vector<double>(Ny) ),
    jx1( Nx+1, std::vector<double>(Ny+1) ),
    jy1( Nx+1, std::vector<double>(Ny+1) )
{

    d = temp/gamma;
    D = temp/Gamma;
    for(int xi=0; xi<Nx; ++xi) {
        x[xi] = 0.5*dx + dx*xi;

        for(int yi=0; yi<Ny; ++yi) {
            y[yi] = 0.5*dy + dy*yi;

            r[xi][yi] =  rInit[xi][yi];
            s[xi][yi] =  sInit[xi][yi];
        }
    }
}


void System::next_time()
{
    next_flux();
    next_prob();
}


void System::next_flux()
{


    // Boundary conditions
    for(int xi=0; xi<Nx+1; ++xi) {
           jx0[xi][0] = 0; 
           jx0[xi][Ny] = 0; 

           jx1[xi][0] = 0; 
           jx1[xi][Ny] = 0; 

           jy0[xi][0] = 0; 
           jy0[xi][Ny] = 0; 

           jy1[xi][0] = 0; 
           jy1[xi][Ny] = 0; 

    }

    for(int yi=0; yi<Ny+1; ++yi) {
           jx0[0][yi] = 0; 
           jx0[Nx][yi] = 0; 

           jx1[0][yi] = 0; 
           jx1[Nx][yi] = 0; 

           jy0[0][yi] = 0; 
           jy0[Nx][yi] = 0;

           jy1[0][yi] = 0; 
           jy1[Nx][yi] = 0; 
    }

    for(int xi=0; xi<Nx-1; ++xi) {
        double xmid = x[xi] + 0.5*dx;
        for(int yi=0; yi<Ny-1; ++yi) {
            double ymid = y[yi] + 0.5*dy;
            jx0[xi][yi]  = 0;
            //jx0[xi][yi] -= swimspeed.v(xmid)*(s[xi][yi] + s[xi+1][yi])/2; 
            //jx0[xi][yi] -= potential.F(xmid, ymid)*(r[xi][yi] + r[xi+1][yi])/(2*gamma); 
            jx0[xi][yi] -= d*(r[xi+1][yi] - r[xi][yi])/dx; 

            jy0[xi][yi]  = 0;
            //jy0[xi][yi] -= potential.F(xmid, ymid)*(r[xi][yi] + r[xi][yi+1])/(2*Gamma);
            jy0[xi][yi] -= D*(r[xi][yi+1] - r[xi][yi])/dy; 

        }
    }


}

void System::next_prob()
{


    for(int xi=0; xi<Nx; ++xi) {
        for(int yi=0; yi<Ny; ++yi) {
            r[xi][yi] -= dy*dt*( jx0[xi+1][yi] - jx0[xi][yi] )/dx;
            r[xi][yi] -= dx*dt*( jy0[xi][yi+1] - jy0[xi][yi] )/dy;

        }
    }

}


void System::save_r( std::ofstream& out ) const
{
    char sep1 = ';';
    char sep2 = '\n';

    for(int xi=0; xi<Nx; ++xi) {
        for(int yi=0; yi<Ny; ++yi) {
            out << r[xi][yi];
            if(yi < (Ny - 1) ) out << sep1;
        }
        if( (xi < Nx - 1) ) out << sep2;
    }
   
}






 
