
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


    // integration time 
    double dt   = config.read<double>("dt");
    

	double epsilon = config.read<double>("epsilon");
	int Nss = config.read<int>("Nss");


    vector<vector<double> > rInit(Nx, vector<double>(Ny,0));
    //rInit[int(Nx/2.)][int(Ny/2.)] = 1./(dx*dy);
	for(int xi=0;xi<Nx; ++xi) {
		rInit[xi][int(Ny/2.)] = 1./(Lx*dy);
	}


    vector<vector<double> > sInit(Nx, vector<double>(Ny,0));


    Swimspeed swimspeed(v0, vp, vx0, alpha, Lx);
    Potential potential(k);
    System system(swimspeed, potential, gamma, Gamma, temp,
                 Nx, dx, Ny, dy, dt, rInit, sInit);

    double t = 0;
    int ti = 0;

	bool steadyState = false;

	string dirname = "data/";

	
	//vector<double> aList {1/10., 1/6., 2/6., 3/6., 4/6., 5/6.,
	//					1., 6/5., 6/4., 6/3., 6/2., 6., 10.};

	//vector<double> aList {4., 2., 1., 0.5, 0.25};
	//vector<double> aList {.75,1.,1.25,1.5,1.75,2.};
	vector<double> aList {2.,1.75, 1.5, 1.25, 1., 0.75};


	system.save_xy(dirname, "");

	for( unsigned int ai=0;ai<aList.size(); ++ai) {
		t = 0;
		ti = 0;
		system.set_Gamma(aList[ai]*gamma);

		cout << "a=" << aList[ai] << endl;
		cout << "Neumann criterion: \n";	
		cout << "x: \t" << temp*dt/(gamma*dx*dx) << "\n";	
		cout << "y: \t" << temp*dt/(aList[ai]*gamma*dy*dy) << "\n";	
		cout << endl;
	  



		steadyState = false;
		rInit = system.get_r();
		sInit = system.get_s();

		while( steadyState == false) {
			
			system.next_time();
			t += dt;
			ti += 1;

			if( (ti % Nss ) == 0 ) {

				steadyState = steady_state(rInit, sInit, system, epsilon);
				cout << "\tEps = " << epsilon << endl;
				cout << "\tEps = " << steady_state(rInit, sInit, system) << endl;

				rInit = system.get_r();
				sInit = system.get_s();
			}

		}
		cout << "time= " << t << endl;
		cout << "ti = " << ti << endl;
		cout << "------------------------------" <<endl << endl;

		stringstream name;
		name << ai;
		system.save_state(dirname,name.str());
		system.save_state_1d(dirname,name.str());
		system.save_flux(dirname, name.str() );
		
	}


    return 0;
}

