// global.h - global variables
# ifndef GLOBAL_H
# define GLOBAL_H
# endif
#include<iostream>
#include<fstream>
#include<vector>
#include<ctime>
#include <iomanip>
#include<cmath>
#include<cstdlib>
#include<algorithm>
#include<utility>
#include<chrono>
#include"omp.h"
//#include"graphs.h"
using namespace std;


// position, velocity, force, image vectors
extern double* RVF;
//extern double* type;

// void-flag vectors
//extern vector<bool> isvoid, isempty;
//extern vector< vector<int> > neighvoidcell;
//extern vector< vector<bool> > voidadjmat;
extern vector<bool> isbond, isoccupy;
// input LAMMPS configuration-dump file
extern const char *infname;

// global parameters/variables
extern int d; // spatial dimension
extern int ncell[3]; // # of cells along each dimension
extern int nbs;  // nbs = 3 or 5 for 3^d or 5^d neighboring cells
extern int N;	// number of spheres
extern int nthreads; // # of OMP threads
extern int nsteps, nthermo, step;  // total # of MD steps, steps between thermo evals, current step
extern int totcells;  // total # of linked cells
extern double phi;   // volume fraction
extern double D[3], cellsize[3], lo[3],scaledlo[3];   // box dimensions and edges
extern double skin;  // for neighbor lists
extern double rho, vol;  // number density, volume
extern double avgneigh, totavgneigh;
extern double hinv[3][3];
extern int type[10000];
extern double h[3][3];


