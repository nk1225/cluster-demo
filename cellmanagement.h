#include"global.h"
#ifndef CELLMANAGEMENT_H
#define CELLMANAGEMENT_H
#endif

// box classes: create variants for other d
class box3D {
public:
    box3D* neighbor[3][3][3];    // pointers to all neighboring boxes
    vector<int> atomlist;         // elements = atomids of atoms in box
    int pos;    // set by binning routine, is population of box = length of atomlist
};



// add public/private variables/methods for other d
class cell {
   public:
// 3D
      vector< vector< vector<box3D> > > supercell3D;
      vector< vector< vector<int> > > boxpops3D;
      void setupsupercell3D();
      void buildvoidlist3D();

   private:
// 3D
    void boxptrsetup3D(int i, int j, int k, int l, int m, int n);

};


