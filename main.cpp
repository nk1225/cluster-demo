    // IDENTIFY EMPTY CUBIC CELLS AND THEIR CONNECTIVITY
#include"cellmanagement.h"
#include"IO.h"
#include"graphs.h"
#include"sched.h"
#define PI 3.1415926536
// global variables
double* RVF; // here only store positions
//int* type;
int d, N, ncell[3], nthreads, nframes, nbs;
double M;
double cellvol;
int totcells, numempty, numcluster, numdistcluster, maxclustersize;
//, numvoids, numdistvoids, maxvoidsize;
double D[3], cellsize[3], mincellsize, lo[3],scaledlo[3], truecellsize[3];
double hinv[3][3],h[3][3];
//double scaled[3];
int type[10000];
vector<int> clustercelllabel, distinctclusternum;
//vector<int> voidcelllabel, distinctvoidnum;
//vector<bool> isempty, isvoid;
vector<bool> isbond,isoccupy;
vector< vector<bool> > voidadjmat;
//vector< vector<int> > neighvoidcell;
double eps, epshalf, avgneigh, totavgneigh, avgforceevals, totavgforceevals, gam, gfac;
const char *infname;
cell supercell;
IO readwrite;

bool isadj(int h1, int h2);

int main (int argc, const char *argv[]) {
    d = 3;
    N = atoi(argv[1]);
    nframes = atoi(argv[2]);
    mincellsize = atof(argv[3]);
 //   nthreads = atoi(argv[4]);
    infname = argv[4];
     

// various loop and helper variables
    int f, i, j, k, l, m;
   // vector<int> voidcelllabel;
    double avgneighvoids, cellvol, avgdistvoidvol;
    double avgneighcluster, avgdistclustervol;
//    vector<int> conn;

    
// for high resolution timing
    chrono::duration<double> elapsed,totelapsed,nltime,totnltime,verltime,totverltime;
    auto start = chrono::high_resolution_clock::now();

    
// initialize position array
   RVF = (double*) _mm_malloc(d*N*sizeof(double),64);

    
    //type array
   // type = (double*) _mm_malloc(N*sizeof(double),64);
    
    
// open input and output files
    readwrite.openfiles();
    ofstream fout1;    fout1.open("cluster_data");
//    ofstream fout2;    fout2.open("voidadjmats");
//ofstream fout3; fout3.open("isorem");    
  fout1 << setw(2)<<"frame"<<setw(8)<<"distcl"<<setw(8)<<"avgSizecl"<<setw(8)<<"maxSizecl"<<endl;  
// loop over frames
    for (f = 0; f < nframes; f++) {
        auto curr = chrono::high_resolution_clock::now();
        elapsed = curr-start;
        cout << "Frame " << f << ", t = " << elapsed.count() <<  endl;
        readwrite.readframe();  // get coords for this frame
        
// set # of cells along each direction + exact cell sizes
        totcells = 1;
        for (i = 0; i < d; i++) {
            ncell[i] = floor(D[i]/mincellsize);
            cellsize[i] = 1.0/ncell[i];//mincellsize/D[i]; // D[i]/ncell[i];
            truecellsize[i] = D[i]/ncell[i];
	    totcells *= ncell[i];
        }
cellvol =truecellsize[0]*truecellsize[1]*truecellsize[2];
// setup cells 
	distinctclusternum.resize(totcells);
        for (i = 0; i < totcells; i++) distinctclusternum[i] = -1;

    
// set up linked cells
        supercell.setupsupercell3D();
// build the cluster lists
        supercell.buildvoidlist3D();
        cout << "cells lists built." << endl;
        
// count the cluster
        numempty = 0;
        for (i = 0; i < totcells; i++) if (isoccupy[i]) numempty++;
        numcluster = 0;
        for (i = 0; i < totcells; i++) if (isbond[i]) numcluster++;
        
        cout<<"cluster counted" << endl;
        
// store the cluster cell labels
        clustercelllabel.resize(numcluster);
        j = 0;
        for (i = 0; i < totcells; i++) if (isbond[i]) {
            clustercelllabel[j] = i;
	//     fout3<<"j "<< j << " voidcellabel "<< i <<endl;
            j++;
        }
        cout << "numbercluster: "<< numcluster<<endl;
        cout << "cluster stored "<< endl;
// build the cluster adjacency graph
        numdistcluster = 0;  // number of DISTINCT (disconnected) voids
        if (numcluster) {
            Graph voidadjmat(numcluster);   // (numvoids);
            for (i = 0; i < numcluster-1; i++) for (j = i+1; j < numcluster; j++) if (isadj(clustercelllabel[i],clustercelllabel[j])) voidadjmat.addEdge(i,j);
            maxclustersize = 0.0;
            numdistcluster = voidadjmat.connectedComponents();
        }
        cout << "cluster adjacency"<<endl;
        

       // cellvol = truecellsize[0]*truecellsize[1]*truecellsize[2];
        if (!numcluster) avgdistclustervol = 0;
        else avgdistclustervol = (numcluster*cellvol)/numdistcluster;
        //cout  << "   number of cluster = " << numdistcluster << "   max cluster size = " << maxclustersize <<  << "   avg. volume per distinct cluster = " << avgdistclustervol << endl;

// need this kludge for numvoids == 1
        if (numcluster == 1|| numdistcluster == 0) {       numcluster = 0;       numdistcluster = 0;   maxclustersize = 0;    avgdistclustervol = 0;  }
// output void-cell fraction
        fout1 <<setw(2)<< f << setw(8) << numdistcluster  << setw(8) << avgdistclustervol << setw(8) << maxclustersize*cellvol << endl;
        
//        delete[] voidadjmat;
//readwrite.writevoidcoords();
        
    } // end of the f-loop
    
    readwrite.closefiles();
    fout1.close();
//    fout2.close();

// high resolution timing and diagnostics
   auto finish = chrono::high_resolution_clock::now();
   totelapsed = finish-start;
   cout << endl << "Execution took " << totelapsed.count() << " seconds." << endl;
   cout << "Runtime per cell = " << totelapsed.count()/totcells << endl;
   for (i = 0; i < 3; i++) cout << endl;
    
   return 0;
}

// determine whether the void cells stored in voidcelllabel[] are adjacent
bool isadj(int vcl1, int vcl2) {  // vcl1, vcl1 are the void cell labels
    int h1, h2, i1, i2, j1, j2, k1, k2;  // positions of the void cells along each dimension
    h1 = vcl1;                              h2 = vcl2;
// find h1 and h2's z-positions
    k1 = h1%ncell[2];                       k2 = h2%ncell[2];
// find h1 and h2's y-positions
    h1 /= ncell[2];                         h2 /= ncell[2];
    j1 = h1%ncell[0];                       j2 = h2%ncell[0];
// find h1 and h2's x-positions
    h1 /= ncell[0];                         h2 /= ncell[0];
    i1 = h1;                                i2 = h2;
    
//    cout << "Void cells " << vcl1 << " and " << vcl2 << " respectively have i1 j1 k1 = " << i1 << " " << j1 << " " << k1 << " and i2 j2 k2 = " << i2 << " " << j2 << " " << k2 << endl;
    
    int d1, d2, d3;   // cell position differences
    d1 = i2 - i1;
    if (d1 < -ncell[0]/2) d1 -= ncell[0];       else if (d1 > ncell[0]/2) d1 += ncell[0];
    d2 = j2 - j1;
    if (d2 < -ncell[1]/2) d2 -= ncell[1];       else if (d2 > ncell[1]/2) d2 += ncell[1];
    d3 = k2 - k1;
    if (d3 < -ncell[2]/2) d3 -= ncell[2];       else if (d3 > ncell[2]/2) d3 += ncell[2];
// cells are adjacent if d1, d2, d3 are each -1, 0, or 1
    // 4 is the old version 3 is the new.
    if (d1*d1 + d2*d2 + d3*d3 < 3) {
//        cout << "Void cells " << vcl1 << " and " << vcl2 << " respectively have i1 j1 k1 = " << i1 << " " << j1 << " " << k1 << ", i2 j2 k2 = " << i2 << " " << j2 << " " << k2 << ", dx dy dz = " << d1 << " " << d2 << " " << d3 << " , (delta d)^2 = " << d1*d1 + d2*d2 + d3*d3 << " and are adjacent." << endl;
        return true;
    }
    else return false;
}


/*
        if (numvoids > 1) {
            voidadjmat.resize(numvoids);
            for (i = 0; i < numvoids; i++) voidadjmat[i].resize(numvoids);
            for (i = 0; i < numvoids; i++) for (j = 0; j < numvoids; j++) voidadjmat[i][j] = false;
            for (i = 0; i < numvoids-1; i++) for (j = i+1; j < numvoids; j++) {
                voidadjmat[i][j] = isadj(voidcelllabel[i],voidcelllabel[j]);
                voidadjmat[j][i] = voidadjmat[i][j];
            }
// This seems to work correctly!  Only problem is it's unfeasible for large numvoids
            for (i = 0; i < numvoids; i++) {
                for (j = 0; j < numvoids; j++) fout2 << voidadjmat[i][j] << " ";
                fout2 << endl;
            }
            fout2 << endl;
            avgneighvoids = 0;
            for (i = 0; i < numvoids; i++) for (j = 0; j < numvoids; j++)   avgneighvoids += (1.0*int(voidadjmat[i][j]))/numvoids;
            cout << "Average void has " << avgneighvoids << " neighbors." << endl;
        }
*/
