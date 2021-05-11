/*! \file cc85.cpp
 *  \brief Definitions of functions needed to run CC85 test. 
 *  Functions are members of the Grid3D class. */ 

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

/*! \fn Background(Real n, Real vx, Real vy, Real vz, Real P)
 *  \brief Constant gas properties throughout the box. */
void Grid3D::Background(Real n, Real vx, Real vy, Real vz, Real P)
{
  int i, j, k, id;
  int istart, jstart, kstart, iend, jend, kend;
  Real x_pos, y_pos, z_pos;
  Real mu = 0.6;
  Real rho, T; 

  rho = n*mu*MP / DENSITY_UNIT;
  vx  = vx * 1e5 / VELOCITY_UNIT;
  vy  = vy * 1e5 / VELOCITY_UNIT; 
  vz  = vz * 1e5 / VELOCITY_UNIT;
  P   = P*KB / PRESSURE_UNIT; 

  istart = H.n_ghost;
  iend   = H.nx-H.n_ghost;
  if (H.ny < 1) {
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

    // Now set the initial values of convserved variables
  for(k=kstart; k<kend; k++) {
    for(j=jstart; j<jend; j++) {
      for(i=istart; i<iend; i++) {

        // Get the cell index
        id = i + j*H.nx + k*H.nx*H.ny; 

        // Get the cell-centered position 
        Get_Position(i, j, k, &x_pos, &y_pos, &z_pos);

        // Set the constant initial states 
        C.density[id]     = rho;
        C.momentum_x[id]  = rho*vx; 
        C.momentum_y[id]  = rho*vy; 
        C.momentum_z[id]  = rho*vz; 
        C.Energy[id]      = P/(gama-1.0) + 0.5*rho*(vx*vx + vy*vy + vz*vz);
        #ifdef DE
        C.GasEnergy[id]   = P/(gama-1.0);
        #endif 
      }
    }
  }
}

/*! \fn Inject()
 *  \brief Add energy and mass to the grid */ 
void Grid3D::Inject()
{
  int i, j, k, id;
  Real x_pos, y_pos, z_pos;
  Real r, center_x, center_y, center_z;
  Real E_dot, M_dot, q_dot, Q_dot, E_tot, M_tot;
  Real alpha, beta, SFR, R, V;  
  Real vx, vy, vz; 
  
  center_x = 0.5;
  center_y = 0.5;
  center_z = 0.5;

  alpha    = 0.1;
  beta     = 0.3;
  R        = 0.3; 
  SFR      = 10;

  int n_cells = 0;
   
  V = (4/3)*PI*(R*R*R);
  

  // Start adding energy and mass after 100 kyr 
  if (H.t > 100) {

    // Set bool so dt is recalculated after energy injection
    H.SN = true;
      
    // Define Mdot and Edot in CGS Units. 
    // Mdot [Msun/yr] is first converted to CGS units [g/s] to be converted to code units  
    // Edot is already in CGS units [erg/s] to be converted to code units 
    // Define R in code units (kpc)

    Real Msun, yr;
    Msun = MASS_UNIT;
    yr = TIME_UNIT / 1e3; 

    M_dot = (beta  * SFR  * (Msun / yr))   / (MASS_UNIT/TIME_UNIT);   
    E_dot = (alpha * SFR  *  3.1e41    )   / (MASS_UNIT*VELOCITY_UNIT*VELOCITY_UNIT/TIME_UNIT);   
    
    q_dot = M_dot / V; 
    Q_dot = E_dot / V; 

    vx = 0;
    vy = 0;
    vz = 0;
      
    for (k=H.n_ghost; k<H.nz-H.n_ghost; k++) {
      for (j=H.n_ghost; j<H.ny-H.n_ghost; j++) {
        for (i=H.n_ghost; i<H.nx-H.n_ghost; i++) {

          id = i + j*H.nx + k*H.nx*H.ny; 

          // Get the Centered x, y, and z positions
        
          Get_Position(i, j, k, &x_pos, &y_pos, &z_pos);

          r = sqrt( (x_pos-center_x)*(x_pos-center_x) + (y_pos-center_y)*(y_pos-center_y) + (z_pos-center_z)*(z_pos-center_z) );
          if ( r < R ) {
            // Add the mass and energy 
            C.density[id]    += q_dot * H.dt;            
            C.Energy[id]     += Q_dot * H.dt;
            #ifdef DE
            C.GasEnergy[id]  += Q_dot * H.dt; 
            #endif
            n_cells ++;          
          }
        }
      }
    }
    // Confirm how much mass & energy added
    M_tot = n_cells*H.dx*H.dy*H.dz*q_dot*H.dt;
    E_tot = n_cells*H.dx*H.dy*H.dz*Q_dot*H.dt;
    printf("Mass added: %e\n", M_tot);
    printf("Energy added: %e\n", E_tot);
    printf("volume: %e\n", V);
    printf("R: %e\n", R);
    // printf("qDot: %e\n", q_dot);
    // printf("QDot: %e\n", Q_dot);
  }
}  
