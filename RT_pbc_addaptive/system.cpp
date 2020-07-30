#include "system.h"

#include <fstream>
#include <iostream>
#include <cmath>


System::System( Swimspeed swimspeed, Potential potential,
                double gamma, double Gamma, double temp,
                int Nx, double dx, int Ny, double dy,
                const std::vector<std::vector<double> >& rInit,
				const std::vector<std::vector<double> >& sInit)
: swimspeed(swimspeed), potential(potential),
    gamma(gamma), Gamma(Gamma), temp(temp),
    Nx(Nx), dx(dx), Ny(Ny), dy(dy),
    x(Nx), y(Ny),
    r( Nx, std::vector<double>(Ny) ),
    jx0( Nx+1, std::vector<double>(Ny) ),
    jy0( Nx, std::vector<double>(Ny+1) ),
    s( Nx, std::vector<double>(Ny) ),
    jx1( Nx+1, std::vector<double>(Ny) ),
    jy1( Nx, std::vector<double>(Ny+1) ),
    rD( Nx, std::vector<double>(Ny) ),
    jx0D( Nx+1, std::vector<double>(Ny) ),
    jy0D( Nx, std::vector<double>(Ny+1) ),
    sD( Nx, std::vector<double>(Ny) ),
    jx1D( Nx+1, std::vector<double>(Ny) ),
    jy1D( Nx, std::vector<double>(Ny+1) )

{
    q = Gamma/gamma;
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


void System::next_time(double dt)
{
	next_flux();
	next_prob(dt);
}

double System::next_time_check_error(double dt)
{
    //next_flux_dummy1();
    //next_prob_dummy1(dt/3.);
    //
    //next_flux_dummy2();
    //next_prob_dummy2(dt/3.);

    //next_flux_dummy2();
    //next_prob_dummy2(dt/3.);


    next_flux();
    next_prob(dt);

    return this->error()*dt;

}

void System::next_prob(double dt)
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

    // bulk	
    for(int xi=0; xi<Nx; ++xi) {
        double xmid = x[xi] - 0.5*dx;

        for(int yi=0; yi<Ny; ++yi) {
            double ymid = y[yi] - 0.5*dy;


            // x index to the left
            double xl = ( xi == 0 ) ? Nx - 1 : xi-1;

            // postition of the active particle
            double X = xmid + q*y[yi]/(1+q);


            jx0[xi][yi]  = swimspeed.v(X) * ( s[xl][yi] + s[xi][yi] ) / 2;
            jx0[xi][yi] -= d * ( r[xi][yi] - r[xl][yi] ) / dx;
            jx0[xi][yi] /= 1 + q;

            jx1[xi][yi]  = swimspeed.v(X) * ( r[xl][yi] + r[xi][yi] ) / 2;
            jx1[xi][yi] -= d * ( s[xi][yi] - s[xl][yi] ) / dx;
            jx1[xi][yi] /= 1 + q;

            if( yi != 0) {

                // postition of the active particle
                double X = x[xi] + q*ymid/(1+q);

                jy0[xi][yi]  = swimspeed.v(X) * ( s[xi][yi-1] + s[xi][yi] ) / 2;
                jy0[xi][yi] += (1+1/q) * potential.F(ymid) * ( r[xi][yi-1] + r[xi][yi] ) /(2*gamma);
                jy0[xi][yi] -= (1+1/q)*d*( r[xi][yi] - r[xi][yi-1] )/dy;

                jy1[xi][yi]  = swimspeed.v(X) * ( r[xi][yi-1] + r[xi][yi] ) / 2;
                jy1[xi][yi] += (1+1/q) * potential.F(ymid) * ( s[xi][yi-1] + s[xi][yi] ) /(2*gamma);
                jy1[xi][yi] -= (1+1/q)*d*( s[xi][yi] - s[xi][yi-1] )/dy;

            }
        }
    }

    // reflecting bc for y
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


void System::next_prob_dummy1(double dt)
{
    for(int xi=0; xi<Nx; ++xi) {
        for(int yi=0; yi<Ny; ++yi) {
            rD[xi][yi] = r[xi][yi];
            rD[xi][yi] -= dt*( jx0D[xi+1][yi] - jx0D[xi][yi] )/dx;
            rD[xi][yi] -= dt*( jy0D[xi][yi+1] - jy0D[xi][yi] )/dy;

            sD[xi][yi] = (1 - dt*swimspeed.get_alpha() )*s[xi][yi];
            sD[xi][yi] -= dt*( jx1D[xi+1][yi] - jx1D[xi][yi] )/dx;
            sD[xi][yi] -= dt*( jy1D[xi][yi+1] - jy1D[xi][yi] )/dy;

        }
    }
}

void System::next_prob_dummy2(double dt)
{
    for(int xi=0; xi<Nx; ++xi) {
        for(int yi=0; yi<Ny; ++yi) {
            rD[xi][yi] -= dt*( jx0D[xi+1][yi] - jx0D[xi][yi] )/dx;
            rD[xi][yi] -= dt*( jy0D[xi][yi+1] - jy0D[xi][yi] )/dy;

            sD[xi][yi] -= dt*swimspeed.get_alpha()*sD[xi][yi];
            sD[xi][yi] -= dt*( jx1D[xi+1][yi] - jx1D[xi][yi] )/dx;
            sD[xi][yi] -= dt*( jy1D[xi][yi+1] - jy1D[xi][yi] )/dy;

        }
    }
}

void System::next_flux_dummy1()
{

    // bulk	
    for(int xi=0; xi<Nx; ++xi) {
        double xmid = x[xi] - 0.5*dx;

        for(int yi=0; yi<Ny; ++yi) {
            double ymid = y[yi] - 0.5*dy;


            // x index to the left
            double xl = ( xi == 0 ) ? Nx - 1 : xi-1;

            // postition of the active particle
            double X = xmid + q*y[yi]/(1+q);


            jx0D[xi][yi]  = swimspeed.v(X) * ( s[xl][yi] + s[xi][yi] ) / 2;
            jx0D[xi][yi] -= d * ( r[xi][yi] - r[xl][yi] ) / dx;
            jx0D[xi][yi] /= 1 + q;

            jx1D[xi][yi]  = swimspeed.v(X) * ( r[xl][yi] + r[xi][yi] ) / 2;
            jx1D[xi][yi] -= d * ( s[xi][yi] - s[xl][yi] ) / dx;
            jx1D[xi][yi] /= 1 + q;

            if( yi != 0) {

                // postition of the active particle
                double X = x[xi] + q*ymid/(1+q);

                jy0D[xi][yi]  = swimspeed.v(X) * ( s[xi][yi-1] + s[xi][yi] ) / 2;
                jy0D[xi][yi] += (1+1/q) * potential.F(ymid) * ( r[xi][yi-1] + r[xi][yi] ) /(2*gamma);
                jy0D[xi][yi] -= (1+1/q)*d*( r[xi][yi] - r[xi][yi-1] )/dy;

                jy1D[xi][yi]  = swimspeed.v(X) * ( r[xi][yi-1] + r[xi][yi] ) / 2;
                jy1D[xi][yi] += (1+1/q) * potential.F(ymid) * ( s[xi][yi-1] + s[xi][yi] ) /(2*gamma);
                jy1D[xi][yi] -= (1+1/q)*d*( s[xi][yi] - s[xi][yi-1] )/dy;

            }
        }
    }

    // reflecting bc for y
    for(int xi=0; xi<Nx; ++xi) {
        jy0D[xi][0] = 0;
        jy0D[xi][Ny] = 0;

        jy1D[xi][0] = 0;
        jy1D[xi][Ny] = 0;
    }

	// periodic bc in x
	for(int yi=0; yi<Ny; ++yi) {
		jx0D[Nx][yi] = jx0D[0][yi];
		jx1D[Nx][yi] = jx1D[0][yi];
	}
   


}

void System::next_flux_dummy2()
{

    // bulk	
    for(int xi=0; xi<Nx; ++xi) {
        double xmid = x[xi] - 0.5*dx;

        for(int yi=0; yi<Ny; ++yi) {
            double ymid = y[yi] - 0.5*dy;


            // x index to the left
            double xl = ( xi == 0 ) ? Nx - 1 : xi-1;

            // postition of the active particle
            double X = xmid + q*y[yi]/(1+q);


            jx0D[xi][yi]  = swimspeed.v(X) * ( sD[xl][yi] + sD[xi][yi] ) / 2;
            jx0D[xi][yi] -= d * ( rD[xi][yi] - rD[xl][yi] ) / dx;
            jx0D[xi][yi] /= 1 + q;

            jx1D[xi][yi]  = swimspeed.v(X) * ( rD[xl][yi] + rD[xi][yi] ) / 2;
            jx1D[xi][yi] -= d * ( sD[xi][yi] - sD[xl][yi] ) / dx;
            jx1D[xi][yi] /= 1 + q;

            if( yi != 0) {

                // postition of the active particle
                double X = x[xi] + q*ymid/(1+q);

                jy0D[xi][yi]  = swimspeed.v(X) * ( sD[xi][yi-1] + sD[xi][yi] ) / 2;
                jy0D[xi][yi] += (1+1/q) * potential.F(ymid) * ( rD[xi][yi-1] + rD[xi][yi] ) /(2*gamma);
                jy0D[xi][yi] -= (1+1/q)*d*( rD[xi][yi] - rD[xi][yi-1] )/dy;

                jy1D[xi][yi]  = swimspeed.v(X) * ( rD[xi][yi-1] + rD[xi][yi] ) / 2;
                jy1D[xi][yi] += (1+1/q) * potential.F(ymid) * ( sD[xi][yi-1] + sD[xi][yi] ) /(2*gamma);
                jy1D[xi][yi] -= (1+1/q)*d*( sD[xi][yi] - sD[xi][yi-1] )/dy;

            }
        }
    }

    // reflecting bc for y
    for(int xi=0; xi<Nx; ++xi) {
        jy0D[xi][0] = 0;
        jy0D[xi][Ny] = 0;

        jy1D[xi][0] = 0;
        jy1D[xi][Ny] = 0;
    }

	// periodic bc in x
	for(int yi=0; yi<Ny; ++yi) {
		jx0D[Nx][yi] = jx0D[0][yi];
		jx1D[Nx][yi] = jx1D[0][yi];
	}
   


}

double System::error()
{

    double emaxX = 0;
    double emaxY = 0;
    double e; 
    for(int xi=0;xi<Nx;++xi){
        for(int yi=0;yi<Ny;++yi){
            if( r[xi][yi] > 0) {
                e = fabs(jx0[xi][yi])/r[xi][yi]; 
                if(e > emaxX ) emaxX = e;
                e = fabs(jx1[xi][yi])/r[xi][yi]; 
                if(e > emaxX ) emaxX = e;
                e = fabs(jy0[xi][yi])/r[xi][yi]; 
                if(e > emaxY ) emaxY = e;
                e = fabs(jy1[xi][yi])/r[xi][yi]; 
                if(e > emaxY ) emaxY = e;

            }
        }
    }
    if( emaxX/dx > emaxY/dy) {
        return emaxX/dx;
    } else {
        return emaxY/dy;
    }
    //double e = 0;
    //double de;
    //for(int xi=0;xi<Nx;++xi){
    //    for(int yi=0;yi<Ny;++yi) {
    //        de = r[xi][yi] - rD[xi][yi];
    //        e += de*de;
    //        //de = s[xi][yi] - sD[xi][yi];
    //        //e += de*de;
    //    }
    //}
    //return e/(Nx*Ny);
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

void System::set_init()
{
    for(int xi=0;xi<Nx; ++xi) {
        for(int yi=0;yi<Ny; ++yi) {
            r[xi][yi] = 0;
            s[xi][yi] = 0;
        }
    }
    double q = Gamma/gamma;
    int ymid = Ny/2;
    for(int xi=0;xi<Nx; ++xi) {
        double pe = this->swimspeed.v(x[xi]);
        pe = pe*pe/( this->swimspeed.get_alpha() * temp );
        pe *= gamma/(1+q); 
        double a = potential.F(-1)/( swimspeed.get_alpha() * gamma);
        double e = 1 - q*q/( q + (1+q)*a);
        //e = 1-q;
        r[xi][ymid] = std::pow(1+pe, -e/2. ); 
    }
    double norm = 0;
    for(int xi=0;xi<Nx; ++xi) norm += r[xi][ymid]*dx;
    for(int xi=0;xi<Nx; ++xi) r[xi][ymid]/=norm;
    

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
