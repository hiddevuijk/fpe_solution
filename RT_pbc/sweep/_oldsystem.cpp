#include "system.h"


System::System( Swimspeed swimspeed, Potential potential,
                double gamma, double Gamma, double temp,
                int Nx, double dx, int Ny, double dy, double dt,
                const std::vector<std::vector<double> >& rInit,
				const std::vector<std::vector<double> >& sInit)
: swimspeed(swimspeed), potential(potential),
    gamma(gamma), Gamma(Gamma), temp(temp),
    Nx(Nx), dx(dx), Ny(Ny), dy(dy),
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
	next_flux();
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

    // Boundary conditions y (separation)
    for(int xi=0; xi<Nx; ++xi) {
		jy0[xi][0] = 0;
		jy0[xi][Ny] = 0;

		jy1[xi][0] = 0;
		jy1[xi][Ny] = 0;
    }
   
	// Boundary conditions x (center of mass)
    //for(int yi=0; yi<Ny; ++yi) {
    //    double xmid = x[0]-0.5*dx;

	//	jx0[0][yi]  = swimspeed.v(xmid,y[yi])*(s[Nx-1][yi] + s[0][yi])/2;
	//	jx0[0][yi] += (1/gamma - 1/Gamma)*potential.F(y[yi])*( r[Nx-1][yi] + r[0][yi] )/2;
	//	jx0[0][yi] -= (d+D)*(r[0][yi] - r[Nx-1][yi])/(2*dx);

	//	jx1[0][yi]  = swimspeed.v(xmid, y[yi])*(r[Nx-1][yi] + r[0][yi])/2;
	//	jx1[0][yi] += (1/gamma - 1/Gamma)*potential.F(y[yi])*( s[Nx-1][yi] + s[0][yi] )/2;
	//	jx1[0][yi] -= (d+D)*(s[0][yi] - s[Nx-1][yi])/(2*dx);

	//	if( yi == 0 ) {


	//	} else if( yi == Ny-1 ) {


	//	} else  {
	//		// 0.25 bc  1/2 for avg over 2 x pos. and 1/2 for 2 dy distance
	//		// 1./(8dy) bc. 1/2dy for double pos. step, 1/2 for average over x, x-1
	//		// and 1/2 from equation.

	//		jx0[0][yi] -= (d-D)*( r[0][yi+1]   - r[0][yi-1]   )/(8*dy);
	//		jx0[0][yi] -= (d-D)*( r[Nx-1][yi+1] - r[Nx-1][yi-1] )/(8*dy);

	//		jx1[0][yi] -= (d-D)*( s[0][yi+1]   - s[0][yi-1]   )/(8*dy);
	//		jx1[0][yi] -= (d-D)*( s[Nx-1][yi+1] - s[Nx-1][yi-1] )/(8*dy);

	//	}


	//	//jx0[0][yi] = 0;
	//	//jx1[0][yi] = 0;

	//	jx0[Nx][yi] = jx0[0][yi];
	//	jx1[Nx][yi] = jx1[0][yi];

    //}


    // bulk
    for(int xi=0; xi<Nx; ++xi) {
        double xmid = x[xi]-0.5*dx;
        for(int yi=0; yi<Ny; ++yi) {
            double ymid = y[yi]-0.5*dy;

            if( xi != 0) { // not along the y boundary ( xi == Nx already excluded)
				jx0[xi][yi]  = swimspeed.v(xmid,y[yi])*(s[xi-1][yi] + s[xi][yi])/2;
				jx0[xi][yi] += (1/gamma - 1/Gamma)*potential.F(y[yi])*( r[xi-1][yi] + r[xi][yi] )/2;
				jx0[xi][yi] -= (d+D)*(r[xi][yi] - r[xi-1][yi])/(2*dx);

				jx1[xi][yi]  = swimspeed.v(xmid, y[yi])*(r[xi-1][yi] + r[xi][yi])/2;
				jx1[xi][yi] += (1/gamma - 1/Gamma)*potential.F(y[yi])*( s[xi-1][yi] + s[xi][yi] )/2;
				jx1[xi][yi] -= (d+D)*(s[xi][yi] - s[xi-1][yi])/(2*dx);

				if( yi != 0 and yi != Ny-1)  {
					// 0.25 bc  1/2 for avg over 2 x pos. and 1/2 for 2 dy distance
					// 1./(8dy) bc. 1/2dy for double pos. step, 1/2 for average over x, x-1
					// and 1/2 from equation.
					jx0[xi][yi] -= (d-D)*( r[xi][yi+1]   - r[xi][yi-1]   )/(8*dy);
					jx0[xi][yi] -= (d-D)*( r[xi-1][yi+1] - r[xi-1][yi-1] )/(8*dy);

					jx1[xi][yi] -= (d-D)*( s[xi][yi+1]   - s[xi][yi-1]   )/(8*dy);
					jx1[xi][yi] -= (d-D)*( s[xi-1][yi+1] - s[xi-1][yi-1] )/(8*dy);

				}
				// do yi == 0 and yi == Ny-1 later

            }

            if( yi != 0 ) { // not along the x boundary ( yi == Ny already excluded)

				jy0[xi][yi]  = swimspeed.v(x[xi], ymid)*( s[xi][yi-1] + s[xi][yi])/2;
				jy0[xi][yi] += (1/gamma + 1/Gamma)*potential.F(ymid)*( r[xi][yi-1] + r[xi][yi] )/2;
				jy0[xi][yi] -= (d+D)*(r[xi][yi] - r[xi][yi-1])/(2*dy);


				jy1[xi][yi]  = swimspeed.v(x[xi], ymid)*( r[xi][yi-1] + r[xi][yi])/2;	
				jy1[xi][yi] += (1/gamma + 1/Gamma)*potential.F(ymid)*( s[xi][yi-1] + s[xi][yi] )/2;
				jy1[xi][yi] -= (d+D)*(s[xi][yi] - s[xi][yi-1])/(2*dy);

				if(xi == 0) {

					jy0[xi][yi] -= (d-D)*( r[xi+1][yi-1] - r[Nx-1][yi-1] )/(8*dx);
					jy0[xi][yi] -= (d-D)*( r[xi+1][yi]   - r[Nx-1][yi]   )/(8*dx);

					jy1[xi][yi] -= (d-D)*( s[xi+1][yi-1] - s[Nx-1][yi-1] )/(8*dx);
					jy1[xi][yi] -= (d-D)*( s[xi+1][yi]   - s[Nx-1][yi]   )/(8*dx);

				} else if( xi == Nx-1 ) {

					jy0[xi][yi] -= (d-D)*( r[0][yi-1] - r[xi-1][yi-1] )/(8*dx);
					jy0[xi][yi] -= (d-D)*( r[0][yi]   - r[xi-1][yi]   )/(8*dx);

					jy1[xi][yi] -= (d-D)*( s[0][yi-1] - s[xi-1][yi-1] )/(8*dx);
					jy1[xi][yi] -= (d-D)*( s[0][yi]   - s[xi-1][yi]   )/(8*dx);

				} else {
					jy0[xi][yi] -= (d-D)*( r[xi+1][yi-1] - r[xi-1][yi-1] )/(8*dx);
					jy0[xi][yi] -= (d-D)*( r[xi+1][yi]   - r[xi-1][yi]   )/(8*dx);

					jy1[xi][yi] -= (d-D)*( s[xi+1][yi-1] - s[xi-1][yi-1] )/(8*dx);
					jy1[xi][yi] -= (d-D)*( s[xi+1][yi]   - s[xi-1][yi]   )/(8*dx);
				} 
	
            }

        }
    }
	//  jx for yi == 0 
	double q = (d-D)/(d+D);
    for(int xi=0; xi<Nx; ++xi) {
		double xmid = x[xi] - 0.5*dx;

		// first the 1/4 (d-D) d/dr rho derivative ( from Jr)
		jx0[xi][0]  = q * ( jy0[xi][1] + jy0[xi-1][0] ) / 4;
		jx0[xi][0] += (1-q) * swimspeed.v(xmid, y[0]) * ( s[xi][0] + s[xi-1][0] ) / 4;
		jx0[xi][0] += (1-q) * (1/gamma - 1/Gamma) * potential.F(y[0]) * ( r[xi][0] + r[xi-1][0] ) / 4;
		jx0[xi][0] += (1-q*q) * (d+D) * ( r[xi][0] - r[xi-1][0] ) / (4*dx);

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



