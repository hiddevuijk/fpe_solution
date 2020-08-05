
#include "system.h"
#include "swimspeed.h"
#include "potential.h"

#include "ConfigFile.h"

#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>


using namespace std;

typedef vector<vector<double> > matrix;

int main()
{
    ConfigFile config("input.txt");
    // system size
    double Lx = config.read<double>("Lx");
    int Nx = config.read<int>("Nx");
    double dx = Lx/Nx;

    double Ly = config.read<double>("Ly");
    int Ny = config.read<int>("Ny");
    double dy = Ly/Ny;

    // potential
    double k = config.read<double>("k");

    // swim/ tumble
    double v0 = config.read<double>("v0");
    double vp = config.read<double>("vp");
    double vx0 = config.read<double>("vx0");
    double alpha = config.read<double>("alpha");

    // diffusion
    double gamma = config.read<double>("gamma");
    double Gamma = config.read<double>("Gamma");
    double temp = config.read<double>("temp");

    bool periodic = config.read<bool>("periodic");

    // integration time 
    double dt   = config.read<double>("dt");
    double time = config.read<double>("time");
    
    int Nprint = config.read<int>("Nprint");
    int Nsave = config.read<int>("Nsave");

	cout << "Neumann criterion: \n";	
	cout << "x: \t" << temp*dt/(gamma*dx*dx) << "\n";	
	cout << "y: \t" << temp*dt/(Gamma*dy*dy) << "\n";	
  

    vector<vector<double> > rInit(Nx, vector<double>(Ny,0));
    //rInit[int(Nx/2.)][int(Ny/2.)] = 1./(dx*dy);
    rInit[int(Nx/2.)][int(Ny/2.)] = 1./(dx*dy);
	//for(int xi=0;xi<Nx;++xi) {
	//	for(int yi=0;yi<Ny;++yi) {
	//		rInit[xi][yi] = 1./(Lx*Ly);
	//	}
	//}

    vector<vector<double> > sInit(Nx, vector<double>(Ny,0));


    Swimspeed swimspeed(v0, vp, vx0, alpha);
    Potential potential(k,Lx);
    System system(swimspeed, potential, gamma, Gamma, temp,periodic,
                 Nx, dx, Ny, dy, dt, rInit, sInit);

    double t = 0;
    int ti = 0;
    int ni = 0;

	string dirname = "data/";
    while( t < time ) {

        if( (ti % int( time/(dt*Nprint) ) ) == 0) cout << int(time/dt) << '\t' << ti << endl;
        if( (ti % int( time/(dt*Nsave)  ) ) == 0) {
            stringstream sr;
            sr << dirname+"r" << ni << ".dat";
            ofstream out_r(sr.str());
            system.save_r(out_r);

            stringstream ss;
            ss << dirname+"s" << ni << ".dat";
            ofstream out_s(ss.str());
            system.save_s(out_s);

            ni += 1;

        }

        system.next_time();
        t += dt;
        ti += 1;
    }

    ofstream out_r(dirname+"r.dat");
    system.save_r(out_r);

    ofstream out_s(dirname+"s.dat");
    system.save_s(out_s);

    ofstream out_x(dirname+"x.dat");
    system.save_x(out_x);

    ofstream out_y(dirname+"y.dat");
    system.save_y(out_y);

	ofstream out_jx0(dirname+"jx0.dat");
	system.save_jx0(out_jx0);

	ofstream out_jx1(dirname+"jx1.dat");
	system.save_jx1(out_jx1);

	ofstream out_jy0(dirname+"jy0.dat");
	system.save_jy0(out_jy0);

	ofstream out_jy1(dirname+"jy1.dat");
	system.save_jy1(out_jy1);

    return 0;
}

