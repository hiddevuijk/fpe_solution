#include "system.h"


System::System( Swimspeed swimspeed, Potential potential,
                double gamma, double Gamma, double temp, bool periodic,
                int Nx, double dx, int Ny, double dy, double dt,
                std::vector<std::vector<double> > rInit, std::vector<std::vector<double> > sInit)
: swimspeed(swimspeed), potential(potential),
    gamma(gamma), Gamma(Gamma), temp(temp),
    periodic(periodic), Nx(Nx), dx(dx), Ny(Ny), dy(dy),
    dt(dt), x(Nx), y(Ny),
    r( Nx, std::vector<double>(Ny) ),
    jx0( Nx+1, std::vector<double>(Ny) ),
    jy0( Nx, std::vector<double>(Ny+1) ),
    s( Nx, std::vector<double>(Ny) ),
    jx1( Nx+1, std::vector<double>(Ny) ),
    jy1( Nx, std::vector<double>(Ny+1) )
{

    d = temp/gamma;
    D = temp/Gamma;
    for(int xi=0; xi<Nx; ++xi) {
        x[xi] = -Nx*dx/2. + 0.5*dx + dx*xi; 
        for(int yi=0; yi<Ny; ++yi) {
            y[yi] = -Ny*dy/2. + 0.5*dy + dy*yi;

            r[xi][yi] =  rInit[xi][yi];
            s[xi][yi] =  sInit[xi][yi];
        }
    }
}


void System::next_time()
{
    if( periodic ) {
        next_flux_pbc();
    } else {
        next_flux();
    }

	next_prob();
}

void System::next_prob()
{
    for(int xi=0; xi<Nx; ++xi) {
        for(int yi=0; yi<Ny; ++yi) {
            r[xi][yi] -= dt*( jx0[xi+1][yi] - jx0[xi][yi] )/dx;
            r[xi][yi] -= dt*( jy0[xi][yi+1] - jy0[xi][yi] )/dy;

            s[xi][yi] -= dt*swimspeed.get_alpha()*s[xi][yi];
            s[xi][yi] -= dt*( jx1[xi+1][yi] - jx1[xi][yi] )/dx;
            s[xi][yi] -= dt*( jy1[xi][yi+1] - jy1[xi][yi] )/dy;

        }
    }
}


void System::next_flux()
{

    // Boundary conditions
    for(int xi=0; xi<Nx; ++xi) {
            jy0[xi][0] = 0;
            jy0[xi][Ny] = 0;

            jy1[xi][0] = 0;
            jy1[xi][Ny] = 0;
    }
    
    for(int yi=0; yi<Ny; ++yi) {
            jx0[0][yi] = 0;
            jx0[Nx][yi] = 0;

            jx1[0][yi] = 0;
            jx1[Nx][yi] = 0;
    }


    // bulk
    for(int xi=0; xi<Nx; ++xi) {
        double xmid = x[xi]-0.5*dx;
        for(int yi=0; yi<Ny; ++yi) {
            double ymid = y[yi]-0.5*dy;

            if( xi != 0) {
                jx0[xi][yi]  = swimspeed.v(xmid)*(s[xi-1][yi] + s[xi][yi])/2;
                jx0[xi][yi] -= d*(r[xi][yi] - r[xi-1][yi])/dx;
                jx0[xi][yi] += potential.F(xmid, y[yi])*( r[xi-1][yi] + r[xi][yi])/(2*gamma);

                jx1[xi][yi]  = swimspeed.v(xmid)*(r[xi-1][yi] + r[xi][yi])/2;
                jx1[xi][yi] -= d*(s[xi][yi] - s[xi-1][yi])/dx;
                jx1[xi][yi] += potential.F(xmid, y[yi])*( s[xi-1][yi] + s[xi][yi])/(2*gamma);
            }

            if( yi != 0 ) {
                jy0[xi][yi]  = 0;
                jy0[xi][yi] -= D*(r[xi][yi] - r[xi][yi-1])/dy; 
                jy0[xi][yi] -= potential.F(x[xi], ymid)*( r[xi][yi-1] + r[xi][yi])/(2*Gamma);

                jy1[xi][yi]  = 0;
                jy1[xi][yi] -= D*(s[xi][yi] - s[xi][yi-1])/dy; 
                jy1[xi][yi] -= potential.F(x[xi], ymid)*( s[xi][yi-1] + s[xi][yi])/(2*Gamma);
            }

        }
    }

}

void System::next_flux_pbc()
{

	// periodic forces

	double xmid, ymid;

    // Boundary conditions
	ymid= y[0]-0.5*dy;
	for(int xi=0; xi<Nx; ++xi) {
		xmid = x[xi]-0.5*dx;

		jy0[xi][0]  = 0;
		jy0[xi][0] -= D*(r[xi][0] - r[xi][Ny-1])/dy; 
		jy0[xi][0] -= potential.F_pbc(x[xi], ymid)*( r[xi][Ny-1] + r[xi][0])/(2*Gamma);

		jy1[xi][0]  = 0;
		jy1[xi][0] -= D*(s[xi][0] - s[xi][Ny-1])/dy; 
		jy1[xi][0] -= potential.F_pbc(x[xi], ymid)*( s[xi][Ny-1] + s[xi][0])/(2*Gamma);

		jy0[xi][0] = 0;
		jy1[xi][0] = 0;

		jy0[xi][Ny] = jy0[xi][0];
		jy1[xi][Ny] = jy1[xi][0];


    }


   
	xmid = x[0] - 0.5*dx; 
    for(int yi=0; yi<Ny; ++yi) {
		ymid = y[yi] - 0.5*dx;
		jx0[0][yi]  = swimspeed.v(xmid)*(s[Nx-1][yi] + s[0][yi])/2;
		jx0[0][yi] -= d*(r[0][yi] - r[Nx-1][yi])/dx;
		jx0[0][yi] += potential.F_pbc(xmid, y[yi])*( r[Nx-1][yi] + r[0][yi])/(2*gamma);

		jx1[0][yi]  = swimspeed.v(xmid)*(r[Nx-1][yi] + r[0][yi])/2;
		jx1[0][yi] -= d*(s[0][yi] - s[Nx-1][yi])/dx;
		jx1[0][yi] += potential.F_pbc(xmid, y[yi])*( s[Nx-1][yi] + s[0][yi])/(2*gamma);


		jx0[0][yi] = 0;
		jx1[0][yi] = 0;

		jx0[Nx][yi] = jx0[0][yi];
		jx1[Nx][yi] = jx1[0][yi];
    }


    // bulk
    for(int xi=0; xi<Nx; ++xi) {
        xmid = x[xi]-0.5*dx;
        for(int yi=0; yi<Ny; ++yi) {
            ymid = y[yi]-0.5*dy;

            if( xi != 0) {
                jx0[xi][yi]  = swimspeed.v(xmid)*(s[xi-1][yi] + s[xi][yi])/2;
                jx0[xi][yi] -= d*(r[xi][yi] - r[xi-1][yi])/dx;
                jx0[xi][yi] += potential.F_pbc(xmid, y[yi])*( r[xi-1][yi] + r[xi][yi])/(2*gamma);
                //jx0[xi][yi] += potential.F_pbc(x[xi-1], y[yi])*r[xi-1][yi]/(2*gamma);
				//jx0[xi][yi] += potential.F_pbc(x[xi],y[yi])*r[xi][yi]/(2*gamma);

                jx1[xi][yi]  = swimspeed.v(xmid)*(r[xi-1][yi] + r[xi][yi])/2;
                jx1[xi][yi] -= d*(s[xi][yi] - s[xi-1][yi])/dx;
                jx1[xi][yi] += potential.F_pbc(xmid, y[yi])*( s[xi-1][yi] + s[xi][yi])/(2*gamma);
                //jx1[xi][yi] += potential.F_pbc(x[xi-1], y[yi])*s[xi-1][yi]/(2*gamma);
                //jx1[xi][yi] += potential.F_pbc(x[xi], y[yi])*s[xi][yi]/(2*gamma);
            }

            if( yi != 0 ) {
                jy0[xi][yi]  = 0;
                jy0[xi][yi] -= D*(r[xi][yi] - r[xi][yi-1])/dy; 
                jy0[xi][yi] -= potential.F_pbc(x[xi], ymid)*( r[xi][yi-1] + r[xi][yi])/(2*Gamma);
                //jy0[xi][yi] -= potential.F_pbc(x[xi], y[yi-1])*r[xi][yi-1]/(2*Gamma);
                //jy0[xi][yi] -= potential.F_pbc(x[xi], y[yi])*r[xi][yi]/(2*Gamma);

                jy1[xi][yi]  = 0;
                jy1[xi][yi] -= D*(s[xi][yi] - s[xi][yi-1])/dy; 
                jy1[xi][yi] -= potential.F_pbc(x[xi], ymid)*( s[xi][yi-1] + s[xi][yi])/(2*Gamma);
                //jy1[xi][yi] -= potential.F_pbc(x[xi], y[yi-1])*s[xi][yi-1]/(2*Gamma);
                //jy1[xi][yi] -= potential.F_pbc(x[xi], y[yi])*s[xi][yi]/(2*Gamma);
            }

        }
    }

}


/*
	Save the state of the function
*/

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

void System::save_s( std::ofstream& out ) const
{
    char sep1 = ';';
    char sep2 = '\n';

    for(int xi=0; xi<Nx; ++xi) {
        for(int yi=0; yi<Ny; ++yi) {
            out << s[xi][yi];
            if(yi < (Ny - 1) ) out << sep1;
        }
        if( (xi < Nx - 1) ) out << sep2;
    }
   
}


void System::save_x( std::ofstream& out ) const
{
    char sep = '\n';

    for(int xi=0; xi<Nx; ++xi) {
        out << x[xi];
        if(xi < (Nx - 1) ) out << sep;
    }
   
}

void System::save_y( std::ofstream& out ) const
{
    char sep = '\n';

    for(int yi=0; yi<Ny; ++yi) {
        out << y[yi];
        if(yi < (Ny - 1) ) out << sep;
    }
   
}

void System::save_jx0( std::ofstream& out ) const
{
    char sep1 = ';';
    char sep2 = '\n';

    for(int xi=0; xi<Nx+1; ++xi) {
        for(int yi=0; yi<Ny; ++yi) {
            out << jx0[xi][yi];
            if(yi < (Ny - 1) ) out << sep1;
        }
        if( (xi < Nx ) ) out << sep2;
    }
   
}

void System::save_jx1( std::ofstream& out ) const
{
    char sep1 = ';';
    char sep2 = '\n';

    for(int xi=0; xi<Nx+1; ++xi) {
        for(int yi=0; yi<Ny; ++yi) {
            out << jx1[xi][yi];
            if(yi < (Ny - 1) ) out << sep1;
        }
        if( (xi < Nx ) ) out << sep2;
    }
   
}


void System::save_jy0( std::ofstream& out ) const
{
    char sep1 = ';';
    char sep2 = '\n';

    for(int xi=0; xi<Nx; ++xi) {
        for(int yi=0; yi<Ny+1; ++yi) {
            out << jy0[xi][yi];
            if(yi < Ny  ) out << sep1;
        }
        if( (xi < Nx-1 ) ) out << sep2;
    }
   
}


void System::save_jy1( std::ofstream& out ) const
{
    char sep1 = ';';
    char sep2 = '\n';

    for(int xi=0; xi<Nx; ++xi) {
        for(int yi=0; yi<Ny+1; ++yi) {
            out << jy1[xi][yi];
            if(yi < Ny  ) out << sep1;
        }
        if( (xi < Nx-1 ) ) out << sep2;
    }
   
}



