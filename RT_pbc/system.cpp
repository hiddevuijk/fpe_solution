#include "system.h"

#include <fstream>
#include <iostream>


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

   	
    for(int xi=0; xi<Nx; ++xi) {

		int xl = xi - 1;	
		if( xl < 0 ) xl = Nx-1;

		int xr = xi + 1;
		if( xr == Nx ) xr = 0;

        double xmid = x[xi]-0.5*dx;

        for(int yi=1; yi<Ny; ++yi) {
			// do yi == 0 and yi == Ny-1 later

            double ymid = y[yi]-0.5*dy;

			if( yi < Ny-1 ) {
				jx0[xi][yi]  = swimspeed.v( xmid,y[yi] ) * ( s[xl][yi] + s[xi][yi] )/4;
				jx0[xi][yi] += (1/gamma - 1/Gamma)*potential.F(y[yi])*( r[xl][yi] + r[xi][yi] )/4;
				jx0[xi][yi] -= (d+D)*( r[xi][yi] - r[xl][yi] )/(4*dx);
				jx0[xi][yi] -= (d-D)*( r[xi][yi+1]   - r[xi][yi-1]  )/(16*dy);
				jx0[xi][yi] -= (d-D)*( r[xl][yi+1] - r[xl][yi-1] )/(16*dy);

				jx1[xi][yi]  = swimspeed.v(xmid, y[yi])*(r[xl][yi] + r[xi][yi])/4;
				jx1[xi][yi] += (1/gamma - 1/Gamma)*potential.F(y[yi])*( s[xl][yi] + s[xi][yi] )/4;
				jx1[xi][yi] -= (d+D)*(s[xi][yi] - s[xl][yi])/(4*dx);
				jx1[xi][yi] -= (d-D)*( s[xi][yi+1]   - s[xi][yi-1]   )/(16*dy);
				jx1[xi][yi] -= (d-D)*( s[xl][yi+1] - s[xl][yi-1] )/(16*dy);
			}


			jy0[xi][yi]  = swimspeed.v(x[xi], ymid)*( s[xi][yi-1] + s[xi][yi])/4;
			jy0[xi][yi] += (1/gamma + 1/Gamma)*potential.F(ymid)*( r[xi][yi-1] + r[xi][yi] )/4;
			jy0[xi][yi] -= (d+D)*(r[xi][yi] - r[xi][yi-1])/(4*dy);
			jy0[xi][yi] -= (d-D)*( r[xr][yi-1] - r[xl][yi-1] )/(16*dx);
			jy0[xi][yi] -= (d-D)*( r[xr][yi]   - r[xl][yi]   )/(16*dx);

			jy1[xi][yi]  = swimspeed.v(x[xi], ymid)*( r[xi][yi-1] + r[xi][yi])/4;	
			jy1[xi][yi] += (1/gamma + 1/Gamma)*potential.F(ymid)*( s[xi][yi-1] + s[xi][yi] )/4;
			jy1[xi][yi] -= (d+D)*(s[xi][yi] - s[xi][yi-1])/(4*dy);
			jy1[xi][yi] -= (d-D)*( s[xr][yi-1] - s[xl][yi-1] )/(16*dx);
			jy1[xi][yi] -= (d-D)*( s[xr][yi]   - s[xl][yi]   )/(16*dx);

        }
    }


	//  jx for yi=0,Ny-1
	double q = (d-D)/(d+D);
    for(int xi=0; xi<Nx; ++xi) {
		int xl = xi - 1;	
		if( xl < 0 ) xl = Nx-1;


		double xmid = x[xi] - 0.5*dx;

		// first the 1/4 (d-D) d/dr rho derivative ( from Jr)
		jx0[xi][0]  = q * ( jy0[xi][1] + jy0[xl][1] ) / 4;
		jx0[xi][0] += (1-q) * swimspeed.v(xmid, y[0]) * ( s[xi][0] + s[xl][0] ) / 4;
		//jx0[xi][0] += (1-q) * (1/gamma - 1/Gamma) * potential.F(y[0]) * ( r[xi][0] + r[xl][0] ) / 4;
		jx0[xi][0] += (1-q*q) * (d+D) * ( r[xi][0] - r[xl][0] ) / (4*dx);

		jx1[xi][0]  = q * ( jy1[xi][1] + jy1[xl][1] ) / 4;
		jx1[xi][0] += (1-q) * swimspeed.v(xmid, y[0]) * ( r[xi][0] + r[xl][0] ) / 4;
		//jx1[xi][0] += (1-q) * (1/gamma - 1/Gamma) * potential.F(y[0]) * ( s[xi][0] + s[xl][0] ) / 4;
		jx1[xi][0] += (1-q*q) * (d+D) * ( s[xi][0] - s[xl][0] ) / (4*dx);

		jx0[xi][Ny-1]  = q * ( jy0[xi][Ny-1] + jy0[xl][Ny-1] ) / 4;
		jx0[xi][Ny-1] += (1-q) * swimspeed.v(xmid, y[Ny-1]) * ( s[xi][Ny-1] + s[xl][Ny-1] ) / 4;
		//jx0[xi][Ny-1] += (1-q) * (1/gamma - 1/Gamma) * potential.F(y[Ny-1]) * ( r[xi][Ny-1] + r[xl][Ny-1] ) / 4;
		jx0[xi][Ny-1] += (1-q*q) * (d+D) * ( r[xi][Ny-1] - r[xl][Ny-1] ) / (4*dx);

		jx1[xi][Ny-1]  = q * ( jy1[xi][Ny-1] + jy1[xl][Ny-1] ) / 4;
		jx1[xi][Ny-1] += (1-q) * swimspeed.v(xmid, y[Ny-1]) * ( r[xi][Ny-1] + r[xl][Ny-1] ) / 4;
		//jx1[xi][Ny-1] += (1-q) * (1/gamma - 1/Gamma) * potential.F(y[Ny-1]) * ( s[xi][Ny-1] + s[xl][Ny-1] ) / 4;
		jx1[xi][Ny-1] += (1-q*q) * (d+D) * ( s[xi][Ny-1] - s[xl][Ny-1] ) / (4*dx);

	}

	
	// Boundary conditions y (separation)
    for(int xi=0; xi<Nx; ++xi) {
		jy0[xi][0] = 0;
		jy0[xi][Ny] = 0;

		jy1[xi][0] = 0;
		jy1[xi][Ny] = 0;
    }

	// periodic bc in x
	for(int yi=0; yi<Ny; ++yi) {
		jx0[Nx][yi] = jx0[0][yi];
		jx1[Nx][yi] = jx1[0][yi];
	}
   


}

/*
	get rx, ry
*/
std::vector<double> System::get_rx() const
{

	std::vector<double> rx(Nx,0);
	for( int xi=0; xi<Nx; ++xi) {
		for( int yi=0; yi<Ny; ++yi) {
			rx[xi] += r[xi][yi]*dy;
		}
	}
	return rx;	

}

std::vector<double> System::get_ry() const
{

	std::vector<double> ry(Ny,0);
	for( int xi=0; xi<Nx; ++xi) {
		for( int yi=0; yi<Ny; ++yi) {
			ry[yi] += r[xi][yi]*dx;
		}
	}
	return ry;	

}


/*
	Save the state of the function
*/

void System::save_state( std::string path, std::string name ) const
{

	std::ofstream out_r(path+"r"+name+".dat");	
	save_r(out_r);

	std::ofstream out_s(path+"s"+name+".dat");	
	save_s(out_s);

}
void System::save_state_1d( std::string path, std::string name ) const
{

	std::vector<double> rx = get_rx();
	std::vector<double> ry = get_ry();

	std::ofstream out_rx(path+"rx"+name+".dat");
	for(int xi=0;xi<Nx; ++xi) {
		out_rx << rx[xi];
		if( xi < Nx-1) out_rx << ';';
	}
	out_rx.close();

	std::ofstream out_ry(path+"ry"+name+".dat");
	for(int yi=0;yi<Ny;++yi) {
		out_ry << ry[yi];
		if( yi<Ny-1) out_ry << ';';
	}
	out_ry.close();

}
void System::save_flux( std::string path, std::string name) const
{

	std::ofstream out_jx0(path+"jx0"+name+".dat");
	save_jx0(out_jx0);
	out_jx0.close();

	std::ofstream out_jx1(path+"jx1"+name+".dat");
	save_jx1(out_jx1);
	out_jx1.close();

	std::ofstream out_jy0(path+"jy0"+name+".dat");
	save_jy0(out_jy0);
	out_jy0.close();

	std::ofstream out_jy1(path+"jy1"+name+".dat");
	save_jy1(out_jy1);
	out_jy1.close();


}

void System::save_xy( std::string path, std::string name) const 
{
    std::ofstream out_x(path+"x"+name+".dat");
    save_x(out_x);
	out_x.close();

    std::ofstream out_y(path+"y"+name+".dat");
    save_y(out_y);
	out_y.close();

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


/*
	NONMEMBER FUNCTIONS
*/

bool steady_state( const std::vector<std::vector<double> >& r,
				   const std::vector<std::vector<double> >& s,
				   const System& system, double epsilon)
{

	double diff = 0;
	double dx = system.get_dx();
	double dy = system.get_dy();

	for(int xi=0; xi<system.get_Nx(); ++xi) {
		for(int yi=0; yi<system.get_Ny(); ++yi) {
			diff += dx*dy*fabs( r[xi][yi] - system.get_r(xi,yi) );
			diff += dx*dy*fabs( s[xi][yi] - system.get_s(xi,yi) );
		}
	}

	//std::cout << diff << std::endl;		

	return diff < epsilon;

}
