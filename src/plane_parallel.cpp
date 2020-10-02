/*! \file plane_parallel.cpp
 *  \brief Definitions of functions needed to run the plane-parallel mass injection test.
           Functions are members of the Grid3D class. */

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "global.h"
#include "grid3D.h"
#include "mpi_routines.h"
#include "error_handling.h"
#include <stdio.h>
#include <cmath>
#include <iostream>
#include <fstream>

using namespace std;

/*! \fn void Wind(Real n, Real vx, Real vy, Real vz, Real P)
 *  \brief Constant gas properties througout the box. */
void Grid3D::Wind(Real n, Real vx, Real vy, Real vz, Real P)
{
  int i, j, k, id;
  int istart, jstart, kstart, iend, jend, kend;
  Real x_pos, y_pos, z_pos;
  Real mu = 0.6;
  Real rho, T;

  // First we need to convert those input values to code units 
  // Code units are found in the global.h header
  // First convert our input units to cgs, then divide by the relevant code unit
  rho = n*mu*MP / DENSITY_UNIT;
  vx  = vx * 1e5 / VELOCITY_UNIT;
  vy  = vy * 1e5 / VELOCITY_UNIT;
  vz  = vz * 1e5 / VELOCITY_UNIT;
  P   = P*KB / PRESSURE_UNIT;


  istart = H.n_ghost;
  iend   = H.nx-H.n_ghost;
  if (H.ny > 1) {
    jstart = H.n_ghost;
    jend   = H.ny-H.n_ghost;
  }
  else {
    jstart = 0;
    jend   = H.ny;
  }
  if (H.nz > 1) {
    kstart = H.n_ghost;
    kend   = H.nz-H.n_ghost;
  }
  else {
    kstart = 0;
    kend   = H.nz;
  }

  // set initial values of conserved variables
  for(k=kstart; k<kend; k++) {
    for(j=jstart; j<jend; j++) {
      for(i=istart; i<iend; i++) {

        //get cell index
        id = i + j*H.nx + k*H.nx*H.ny;

        // get cell-centered position
        Get_Position(i, j, k, &x_pos, &y_pos, &z_pos);
        
        // set constant initial states
        C.density[id]    = rho;
        C.momentum_x[id] = rho*vx;
        C.momentum_y[id] = rho*vy;
        C.momentum_z[id] = rho*vz;
        C.Energy[id]     = P/(gama-1.0) + 0.5*rho*(vx*vx + vy*vy + vz*vz);
        #ifdef DE
        C.GasEnergy[id]  = P/(gama-1.0);
        #endif

        // comment this out if you like, I just like to make sure
        // the initial values are sensible
        if (i==istart && j==jstart && k==kstart) {
          n = rho*DENSITY_UNIT / (mu*MP);
          T = P*PRESSURE_UNIT / (n*KB);
          printf("Initial n = %e, T = %e\n", n, T);
        }

      }
    }
  }

}



/*! \fn void Wind_Boundary()
 *  \brief Apply a wind boundary to the -x faces */
void Grid3D::Wind_Boundary()
{
  int i, j, k, id;
  Real x_pos, y_pos, z_pos, r;
  Real n, rho, vx, vy, vz, P, T;
  // for now, I just hard-coded these values in - in theory
  // I should rewrite this function to be able to take parameters
  // from the parameter file, but I'm feeling lazy
  Real mu = 0.6;
  n = 0.005;
  vx = 1000;
  vy = 0.0;
  vz = 0.0;
  P  = 3.15e4;

  // as in Wind IC function, convert these values to code units
  rho = n*mu*MP / DENSITY_UNIT;
  vx = vx*1e5 / VELOCITY_UNIT;
  vy = vy*1e5 / VELOCITY_UNIT;
  vz = vz*1e5 / VELOCITY_UNIT;
  P  = P*KB / PRESSURE_UNIT;

  // set exact boundaries on the -x face
  for (k=0; k<H.nz; k++) {
    for (j=0; j<H.ny; j++) {
      for (i=0; i<H.n_ghost; i++) {

        id = i + j*H.nx + k*H.nx*H.ny;
        // get the (centered) x, y, and z positions at (x,y,z)
        Get_Position(i, j, k, &x_pos, &y_pos, &z_pos);
        // set the conserved quantities
        C.density[id]  = rho;
        C.momentum_x[id] = vx*C.density[id];
        C.momentum_y[id] = vy*C.density[id];
        C.momentum_z[id] = vz*C.density[id];
        C.Energy[id]     = P/(gama-1.0) + 0.5*C.density[id]*(vx*vx + vy*vy + vz*vz);
        #ifdef DE
        C.GasEnergy[id]  = P/(gama-1.0);
        #endif //DE
/*
        if (i==0 && j==0 && k==0) {
          n = rho*DENSITY_UNIT / (mu*MP);
          T = P*PRESSURE_UNIT / (n*KB);
          printf("Boundary n: %e,  T: %e\n", n, T);
        }
*/
      }
    }
  }

}



/*! \fn void Add_Mass()
 *  \brief Add mass to the grid */
void Grid3D::Add_Mass()
{
  int i, j, k, id;
  Real x_pos, y_pos, z_pos;
  Real M_dot, V, rho_dot, M_tot;
  int n_cells = 0;

  // start adding mass after 100 kyr
  if (H.t > 100) {

    // define Mdot in code units
    // these are set in global.h, currently Msun / kyr
    M_dot = 10.0;

    // define the region within which you want to deposit mass
    Real xmin, xmax, ymin, ymax, zmin, zmax;
    xmin = 0.5;
    xmax = 0.6;
    ymin = 0.0;
    ymax = 1.0;
    zmin = 0.0;
    zmax = 1.0;

    V = (xmax - xmin)*(ymax - ymin)*(zmax - zmin);

    // convert Mdot to a mass injection rate per volume
    rho_dot = M_dot / V;

    // loop over all the cells
    for (k=H.n_ghost; k<H.nz-H.n_ghost; k++) {
      for (j=H.n_ghost; j<H.ny-H.n_ghost; j++) {
        for (i=H.n_ghost; i<H.nx-H.n_ghost; i++) {

          id = i + j*H.nx + k*H.nx*H.ny;
          // get the (centered) x, y, and z positions at (x,y,z)
          Get_Position(i, j, k, &x_pos, &y_pos, &z_pos);
          // if the position is in the range where you want to add mass, add mass
          if (x_pos > xmin && x_pos < xmax && y_pos > ymin && y_pos < ymax && z_pos > zmin && z_pos < zmax) {

            // add the mass
            C.density[id]  += rho_dot*H.dt;

            n_cells ++;
          }
        }
      }
    }

    // confirm how much mass we added
    M_tot = n_cells*H.dx*H.dy*H.dz*rho_dot*H.dt;
    printf("Mass added: %e\n", M_tot);

  }

}






