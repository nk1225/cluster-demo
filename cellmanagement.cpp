#include "cellmanagement.h"

int pointvec3D[3];

// 3D assignment of pointers to neighboring boxes (respecting PBCs)
// generalized to having different numbers of cells along different directions
void cell::boxptrsetup3D(int i, int j, int k, int l, int m, int n) {
    int o, p, q;
    // PBCs expressed here: (l, m, n) is the vector from the pointer
    if ((i + l) == -1) o = ncell[0]-1;
    else if ((i + l) == ncell[0]) o = 0;                // handles x rollover
    else o = i + l;
    pointvec3D[0] = o;
    if ((j + m) == -1) p = ncell[1]-1;
    else if ((j + m) == ncell[1]) p = 0;                // handles y rollover
    else p = j + m;
    pointvec3D[1] = p;
    if ((k + n) == -1) q = ncell[2]-1;                    // handles z rollover
    else if ((k + n) == ncell[2]) q = 0;
    else q = k + n;
    pointvec3D[2] = q;
   // cout << " boxsetup finish" <<endl;
    return;
}


// set up 3D supercell, generalized to having different numbers of cells along different directions
void cell::setupsupercell3D() {
    int i, j, k, l, m, n;
    supercell3D.clear();
    supercell3D.resize(ncell[0]);    boxpops3D.resize(ncell[0]);
    for (i = 0; i < ncell[0]; i++) { supercell3D[i].resize(ncell[1]);   boxpops3D[i].resize(ncell[1]); }
    for (i = 0; i < ncell[0]; i++) for (j = 0; j < ncell[1]; j++) { supercell3D[i][j].resize(ncell[2]);   boxpops3D[i][j].resize(ncell[2]); }
    
// link cells; each cell has a pointer to each of its neighbors, 3D-pbcs assumed
    for (i = 0; i < ncell[0]; i++) for (j = 0; j < ncell[1]; j++) for (k = 0; k < ncell[2]; k++) {
       for (l = -1; l < 2; l++) for (m = -1; m < 2; m++) for (n = -1; n < 2; n++) {
           boxptrsetup3D(i,j,k,l,m,n);
           supercell3D[i][j][k].neighbor[l+1][m+1][n+1] = &supercell3D[pointvec3D[0]][pointvec3D[1]][pointvec3D[2]];
       }
    }

    isoccupy.resize(totcells);       isbond.resize(totcells);
//    neighvoidcell.resize(totcells);
//    for (i = 0; i < totcells; i++) neighvoidcell[i].resize(26);
    
    
    cout << "Supercell setup completed." << endl;
    return;
}


// build or rebuild the void and void-cell neighbor lists
// should probably define a bool isempty(totcells) so isempty(i) = true if cell i is empty
// but isvoid(i) is true only if cell i INTERSECTS no atom cores
// e.g. unit-cube cells are voids are only if they are empty and no atoms lie within .5 of their boundaries
void cell::buildvoidlist3D() {
    int i, j, k, l, m, n, p, r, u, v,q;
    //int tox,toy,toz;
    double scaled[d],scadeldel1[v],scadeldel2[d];
    int thiscell, thatcell, thisneigh;
    vector <vector <int>> box;
    box.resize(d);
    for (i = 0; i < d; i++) box[i].resize(N);
    int to[d];    // used to place atoms in bins
    double thislo[d], thishi[d];    // opposite corners of the current cell
    double del1[d], del2[d], delsq;   // distances of other atoms to the current cell's corners (see below)
    
// Put atoms in boxes, set boxpops
    for (i = 0; i < ncell[0]; i++) for (j = 0; j < ncell[1]; j++) for (k = 0; k < ncell[2]; k++) boxpops3D[i][j][k] = 0;
    for (u = 0; u < N; u++) { // should switch order of these loops?
        for (v = 0; v < d; v++) {
            // introduce scaled position for triclinic box
            scaled[v] = hinv[v][0]*RVF[d*u] + hinv[v][1]*RVF[d*u+1]+hinv[v][2]*RVF[d*u+2];
            
            to[v] = static_cast <int> ((scaled[v] - scaledlo[v])/cellsize[v]);
           
           if (to[v] < 0) { to[v] = ncell[v]-1; RVF[3*u+v] += D[v]; }
            else if (to[v] >= ncell[v]) { to[v] = 0; RVF[3*u+v] -= D[v]; }
           box[v][u] = to[v];
            
        }
        boxpops3D[to[0]][to[1]][to[2]] += 1;
    }
     // set up atom lists of length boxpops
    for (i = 0; i < ncell[0]; i++) for (j = 0; j < ncell[1]; j++) for (k = 0; k < ncell[2]; k++) {
        supercell3D[i][j][k].atomlist.resize(boxpops3D[i][j][k]);
        supercell3D[i][j][k].pos = 0;
    }
// second step of putting atoms in boxes
    for (u = 0; u < N; u++) {
        supercell3D[box[0][u]][box[1][u]][box[2][u]].atomlist[supercell3D[box[0][u]][box[1][u]][box[2][u]].pos] = u;
        supercell3D[box[0][u]][box[1][u]][box[2][u]].pos++;
    }
    
// set empty-cell flags, initialize void-cell flags
    for (i = 0; i < ncell[0]; i++) for (j = 0; j < ncell[1]; j++) for (k = 0; k < ncell[2]; k++) {
        thiscell = i*ncell[1]*ncell[2] + j*ncell[2] + k;
        if (!boxpops3D[i][j][k]) {  isoccupy[thiscell] = false;       isbond[thiscell] = false;  }
        else {  isoccupy[thiscell] = true;       isbond[thiscell] = true;  }
    }

// build void and neighbor lists using a d-deep loop over the cells
// empty cells are voids if they also INTERSECT no cores of atoms in their neighboring cells
//    if (nthreads > 1) omp_set_num_threads(nthreads);
// parallelizing the outermost loop (over the cells in the first dimension) is simple and effective
    __assume_aligned(RVF,64);
 //   __assume_aligned(type,64);
//    #pragma omp parallel for private(j,k,l,m,n,thiscell,thatcell,thisneigh)
    for (i = 0; i < ncell[0]; i++) for (j = 0; j < ncell[1]; j++) for (k = 0; k < ncell[2]; k++) {
 //       thiscell = i*ncell[0]*ncell[1] + j*ncell[1] + k;
        thiscell = i*ncell[1]*ncell[2] + j*ncell[2] + k;
        if (isoccupy[thiscell]) {
// define opposite corners of the current cell so can test whether atom cores in other cells intersect it
           // changed from triclinic box
            thislo[0] = scaledlo[0] + i*cellsize[0]; thislo[1] = scaledlo[1] + j*cellsize[1];
            thislo[2] = scaledlo[2] + k*cellsize[2];
            thishi[0] = scaledlo[0] + (i+1)*cellsize[0];      thishi[1] = scaledlo[1] + (j+1)*cellsize[1];      thishi[2] = scaledlo[2] + (k+1)*cellsize[2];
            
// loop over cells neighboring supercell[i][j][k]; check whether any atom cores from neighboring cells intersect thiscell
// see https://stackoverflow.com/questions/4578967/cube-sphere-intersection-test
            for (l = 0; l < 3; l++) for (m = 0; m < 3; m++) for (n = 0; n < 3; n++) {
                for (p = 0; p < supercell3D[i][j][k].neighbor[l][m][n]->pos; p++) {
                    r = supercell3D[i][j][k].neighbor[l][m][n]->atomlist[p];  // index of the atom testing for overlap with the current cell
                    
                    // type = 1 radius = .5
                    if(type[r] == 1){
                    delsq = .25;  // assume sphere core radius = .5
                    for (v = 0; v < d; v++) {
                        scaled[v] = hinv[v][0]*RVF[d*r] + hinv[v][1]*RVF[d*r+1]+hinv[v][2]*RVF[d*r+2];
                         del1[v] = scaled[v] - thislo[v];
                        if (del1[v] < -1.0/2.0) del1[v] += 1;       else if (del1[v] > 1.0/2.0) del1[v] -= 1;
                        
                        
                        del2[v] = scaled[v] - thishi[v];
                        if (del2[v] < -1.0/2.0) del2[v] += 1;       else if (del2[v] > 1.0/2.0) del2[v] -= 1;
		    }
                   for(v = 0; v<d; v++){
		     
		      	   scadeldel1[v] = h[v][0]*del1[0] + h[v][1]*del1[1] + h[v][2]*del1[2];
                            scadeldel2[v] = h[v][0]*del2[0] + h[v][1]*del2[1] + h[v][2]*del2[2];
                      
			   
			    

                        if (scadeldel1[v] < 0) delsq -= scadeldel1[v]*scadeldel1[v];
                            else if (scadeldel2[v] > 0) delsq -= scadeldel2[v]*scadeldel2[v];
                           
		   } 
// atom r's core intersects the current cell if delsq > 0
                    if (delsq > 0) {    isoccupy[thiscell] = true;      goto nextcell;   }
                }
                
                // type = 2 radius = .7
                else if(type[r] == 2) {
                    delsq = .49;  // assume sphere core radius = .7
                    for (v = 0; v < d; v++) {
                        scaled[v] = hinv[v][0]*RVF[d*r] + hinv[v][1]*RVF[d*r+1]+hinv[v][2]*RVF[d*r+2];
                        del1[v] = scaled[v] - thislo[v];
                       // del1[v] = RVF[d*r+v] - thislo[v];
                        
                        if (del1[v] < -1.0/2.0) del1[v] += 1;       else if (del1[v] > 1.0/2.0) del1[v] -= 1;
                        del2[v] = scaled[v] - thishi[v];
                       // del2[v] = RVF[d*r+v] - thishi[v];
                        if (del2[v] < -1.0/2.0) del2[v] += 1;       else if (del2[v] > 1.0/2.0) del2[v] -= 1;
                        }
              for(v = 0; v<d; v++){    
	      	    scadeldel1[v] = h[v][0]*del1[0] + h[v][1]*del1[1] + h[v][2]*del1[2];
                        scadeldel2[v] = h[v][0]*del2[0] + h[v][1]*del2[1] + h[v][2]*del2[2];
                    
		//	fout2<<"atom "<<r<<"v "<<v<<"del1 "<<scadeldel1[v]<<endl;
                  //        fout2<<"atom "<<r<<"v "<<v<<"del2 "<<scadeldel2[v]<<endl;
	      
			  if (scadeldel1[v] < 0) delsq -= scadeldel1[v]*scadeldel1[v];
                        else if (scadeldel2[v] > 0) delsq -= scadeldel2[v]*scadeldel2[v];
			
                    
	      }    
		    // atom r's core intersects the current cell if delsq > 0

                    if (delsq > 0) {    isoccupy[thiscell] = true;      goto nextcell;   }
                }
                }
                nextcell: ;
            } // end of the lmn-loop
        } // end of the isempty[thiscell] conditional
    } // end of the ijk-loop
    // deparallelize
//    if (nthreads > 1) omp_set_num_threads(1);

    return;
}
