
#include "system.h"
#include "swimspeed.h"
#include "potential.h"

#include "ConfigFile.h"

#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <filesystem>


using namespace std;


typedef vector<vector<double> > matrix;

int main()
{

	// Read variable values
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

    // integration time 
    double dt   = config.read<double>("dt");
    double time = config.read<double>("time");
    
    //int Nprint = config.read<int>("Nprint");
    double Tsave = config.read<double>("Tsave");

	bool saveXY = config.read<bool>("saveXY");

	double epsilon = config.read<double>("epsilon");


	double Tss = config.read<double>("Tss");

	string dirname = config.read<string>("dirname");
    bool restart = config.read<bool>("restart");
    
	cout << "Neumann criterion: \n";	
    cout << "x: \t" << temp*dt*(1+gamma/Gamma)/(gamma*dx*dx) << "\n";
    cout << "y: \t" << temp*dt/(gamma*dy*dy*(1+Gamma/gamma)) << "\n";
  

    vector<vector<double> > rInit(Nx, vector<double>(Ny,0));
    //rInit[int(Nx/2.)][int(Ny/2.)] = 1./(dx*dy);
	for(int xi=0;xi<Nx; ++xi) {
		rInit[xi][int(Ny/2.)] = 1./(Lx*dy);
	}
	

    vector<vector<double> > sInit(Nx, vector<double>(Ny,0));

    
    Swimspeed swimspeed(v0, vp, vx0, alpha, Lx);
    Potential potential(k);
    System system(swimspeed, potential, gamma, Gamma, temp,
                 Nx, dx, Ny, dy, rInit, sInit);

    if( restart ) {
        system.read_init(dirname);
        dirname = config.read<string>("new_dirname");
        rInit = system.get_r();
        sInit = system.get_s();
    } else {
        system.set_init();
    }

    dirname += "/";

	vector<double> tList;
    double t = 0;
    int ti = 0;
    int ni = 0;

	bool steadyState = false;
    double steadyStateError;


	tList.push_back(0);

	system.save_state(dirname, "0");
	system.save_state_1d(dirname, "0");


   
    if(dx < dy) {
        cout << (v0+vp)*dt/dx << endl;
    } else {
        cout << (v0+vp)*dt/dy << endl;
    }
         
    std::vector<double> dtList;


    double dTimesave = 0;
    double dTss = 0;
    while( t < time  and steadyState == false ) {
        
        system.next_time(dt);
        t += dt;
        dTimesave += dt;
        dTss += dt;
        ti += 1;
        dtList.push_back(dt);

        if( dTimesave > Tsave ) {
            cout <<"----------"<< endl;
            cout << "ni\t" << ni<< endl;
            cout << "t\t" << t << endl;
            cout << "dt\t"<< dt << endl;
            cout << "eps\t" << epsilon << endl;
            cout << "err\t" << steadyStateError << endl;
            cout <<"----------"<< endl;
            dTimesave -= Tsave;

            ni += 1;

			stringstream name;
			name << ni;
			if(saveXY) {
				system.save_state(dirname,name.str());
			}
			system.save_state_1d(dirname,name.str());

			tList.push_back(t);
        }
		if( dTss > Tss ) {
            dTss = 0;
			steadyStateError = steady_state_error(rInit, sInit, system);
            steadyState = steadyStateError < epsilon;

			rInit = system.get_r();
			sInit = system.get_s();
		}
    }
    cout << t << endl << endl;

	system.save_state(dirname, "");
	system.save_state_1d(dirname, "");
	system.save_flux(dirname, "");
	system.save_xy(dirname, "");
	

	// save time
	ofstream out_time(dirname+"t.dat");
	for(unsigned int i=0;i<tList.size(); i++) {
		out_time << tList[i];
		if(i<tList.size()-1) out_time << ';';
	}
	// save dtList
	ofstream out_dt(dirname+"dt.dat");
	for(unsigned int i=0;i<dtList.size(); i++) {
		out_dt << dtList[i];
		if(i<dtList.size()-1) out_dt << ';';
	}


    return 0;
}

