
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
    double Tsave = config.read<double>("Nsave");

	bool saveXY = config.read<bool>("saveXY");

	double epsilon = config.read<double>("epsilon");

    double dxEmax = config.read<double>("dxEmax");
    double NEmax = config.read<double>("NEmax");


	int Nss = config.read<int>("Nss");

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

    system.set_init();

	vector<double> tList;
    double t = 0;
    int ti = 0;
    int ni = 0;

	bool steadyState = false;

	string dirname = "data/";

	tList.push_back(0);

	system.save_state(dirname, "0");
	system.save_state_1d(dirname, "0");


   
    double q = Gamma/gamma; 
    double dtNeumannR = dx*dx*gamma*(1+q)/(2*temp);
    double dtNeumannr = dy*dy*gamma*q/(2*temp*(1+q));
    double dtNeumann =  dtNeumannR < dtNeumannr ? dtNeumannR : dtNeumannr;

    if(dx < dy) {
        cout << (v0+vp)*dt/dx << endl;
    } else {
        cout << (v0+vp)*dt/dy << endl;
    }
         
    dt = dtNeumann*NEmax/10;
    std::vector<double> dtList;
    int ncheckError =10;
    double dtmin = 1.e-4;

    double dTimesave = 0;
    double e = 0;
    while( t < time  and steadyState == false ) {
        
        if( ti % ncheckError == 0  ){
    
            e = system.next_time_check_error(dt);
            t += dt;
            ti += 1;
            dTimesave += dt;
            if( e > dxEmax ) {
                dt /= 1.5;
            } else{
                dt *= 1.020;
            }
            if( dt > dtNeumann*NEmax ) dt = dtNeumann*NEmax;
            if( dt < dtmin) dt = dtmin;
            
        } else {
            system.next_time(dt);
            t += dt;
            dTimesave += dt;
            ti += 1;
        }
        dtList.push_back(dt);

        if( dTimesave > Tsave ) {
            cout <<"----------"<< endl;
            cout << "ni\t" << ni << endl;
            cout << "t\t" << t << endl;
            cout << "dt\t"<< dt << endl;
            cout << "eN\t" << dtNeumann*NEmax << endl;
            cout << "dx\t" << e << endl;
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
		if( (ti % Nss ) == 0 ) {

			steadyState = steady_state(rInit, sInit, system, epsilon);

			rInit = system.get_r();
			sInit = system.get_s();
		}

    }

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

