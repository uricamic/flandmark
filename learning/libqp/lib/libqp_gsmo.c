/*-----------------------------------------------------------------------
 * libqp_gsmo.c: implementation of the Generalized SMO algorithm.
 *
 * DESCRIPTION
 *  The library provides function which solves the following instance of
 *  a convex Quadratic Programming task:
 *
 *  min QP(x) := 0.5*x'*H*x + f'*x  
 *   x                                      
 *
 *   s.t.    a'*x = b 
 *           LB[i] <= x[i] <= UB[i]   for all i=1..n
 *
 * A precision of the found solution is controlled by the input argument
 * TolKKT which defines tightness of the relaxed Karush-Kuhn-Tucker 
 * stopping conditions.
 *
 * INPUT ARGUMENTS
 *  get_col   function which returns pointer to the i-th column of H.
 *  diag_H [double n x 1] vector containing values on the diagonal of H.
 *  f [double n x 1] vector.
 *  a [double n x 1] Vector which must not contain zero entries.
 *  b [double 1 x 1] Scalar.
 *  LB [double n x 1] Lower bound; -inf is allowed.
 *  UB [double n x 1] Upper bound; inf is allowed.
 *  x [double n x 1] solution vector; must be feasible.
 *  n [uint32_t 1 x 1] dimension of H.
 *  MaxIter [uint32_t 1 x 1] max number of iterations.
 *  TolKKT [double 1 x 1] Tightness of KKT stopping conditions.
 *  print_state  print function; if == NULL it is not called.
 *
 * RETURN VALUE
 *  structure [libqp_state_T]
 *   .QP [1x1] Primal objective value.
 *   .exitflag [1 x 1] Indicates which stopping condition was used:
 *     -1  ... not enough memory
 *      0  ... Maximal number of iterations reached: nIter >= MaxIter.
 *      4  ... Relaxed KKT conditions satisfied. 
 *   .nIter [1x1] Number of iterations.
 *
 * REFERENCE
 *  S.S. Keerthi, E.G. Gilbert. Convergence of a generalized SMO algorithm 
 *   for SVM classier design. Technical Report CD-00-01, Control Division, 
 *   Dept. of Mechanical and Production Engineering, National University 
 *   of Singapore, 2000. 
 *   http://citeseer.ist.psu.edu/keerthi00convergence.html  
 *
 *
 * Copyright (C) 2006-2008 Vojtech Franc, xfrancv@cmp.felk.cvut.cz
 * Center for Machine Perception, CTU FEL Prague
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public 
 * License as published by the Free Software Foundation; 
 * Version 3, 29 June 2007
 *-------------------------------------------------------------------- */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <limits.h>

#include "libqp.h"

libqp_state_T libqp_gsmo_solver(const double* (*get_col)(uint32_t),
            double *diag_H,
            double *f,
            double *a,
            double b,
            double *LB,
            double *UB,
            double *x,
            uint32_t n,
            uint32_t MaxIter,
            double TolKKT,
            void (*print_state)(libqp_state_T state))
{
  double *col_u;
  double *col_v;
  double *Nabla;
  double minF_up;
  double maxF_low;
  double tau;
  double F_i;
  double tau_ub, tau_lb;
  double Q_P;
  uint32_t i, j;
  uint32_t u, v;
  libqp_state_T state;

  Nabla = NULL;

  /* ------------------------------------------------------------ */
  /* Initialization                                               */
  /* ------------------------------------------------------------ */

  /* Nabla = H*x + f is gradient*/
  Nabla = (double*)LIBQP_CALLOC(n, sizeof(double));
  if( Nabla == NULL )
  {
	  state.exitflag=-1;
	  goto cleanup;
  }
 
  /* compute gradient */
  for( i=0; i < n; i++ ) 
  {
    Nabla[i] += f[i];
    if( x[i] != 0 ) {
      col_u = (double*)get_col(i);      
      for( j=0; j < n; j++ ) {
          Nabla[j] += col_u[j]*x[i];
      }
    }
  }

  if( print_state != NULL) 
  {
     state.QP = 0;
     for(i = 0; i < n; i++ ) 
        state.QP += 0.5*(x[i]*Nabla[i]+x[i]*f[i]); 

     print_state( state );
  }


  /* ------------------------------------------------------------ */
  /* Main optimization loop                                       */
  /* ------------------------------------------------------------ */

  state.nIter = 0;
  state.exitflag = 100;
  while( state.exitflag == 100 ) 
  {
    state.nIter ++;     

    /* find the most violating pair of variables */
    minF_up = LIBQP_PLUS_INF;
    maxF_low = -LIBQP_PLUS_INF;
    for(i = 0; i < n; i++ ) 
    {
      
      F_i = Nabla[i]/a[i];

      if(LB[i] < x[i] && x[i] < UB[i]) 
      { /* i is from I_0 */
        if( minF_up > F_i) { minF_up = F_i; u = i; }
        if( maxF_low < F_i) { maxF_low = F_i; v = i; }
      } 
      else if((a[i] > 0 && x[i] == LB[i]) || (a[i] < 0 && x[i] == UB[i])) 
      { /* i is from I_1 or I_2 */
        if( minF_up > F_i) { minF_up = F_i; u = i; }
      }
      else if((a[i] > 0 && x[i] == UB[i]) || (a[i] < 0 && x[i] == LB[i])) 
      { /* i is from I_3 or I_4 */
        if( maxF_low < F_i) { maxF_low = F_i; v = i; }
      }
    }

    /* check KKT conditions */
    if( maxF_low - minF_up <= TolKKT )
      state.exitflag = 4;
    else 
    {
      /* SMO update of the most violating pair */
      col_u = (double*)get_col(u);
      col_v = (double*)get_col(v);

      if( a[u] > 0 ) 
         { tau_lb = (LB[u]-x[u])*a[u]; tau_ub = (UB[u]-x[u])*a[u]; }
      else
         { tau_ub = (LB[u]-x[u])*a[u]; tau_lb = (UB[u]-x[u])*a[u]; }

      if( a[v] > 0 )
         { tau_lb = LIBQP_MAX(tau_lb,(x[v]-UB[v])*a[v]); tau_ub = LIBQP_MIN(tau_ub,(x[v]-LB[v])*a[v]); }
      else
         { tau_lb = LIBQP_MAX(tau_lb,(x[v]-LB[v])*a[v]); tau_ub = LIBQP_MIN(tau_ub,(x[v]-UB[v])*a[v]); }

      tau = (Nabla[v]/a[v]-Nabla[u]/a[u])/
            (diag_H[u]/(a[u]*a[u]) + diag_H[v]/(a[v]*a[v]) - 2*col_u[v]/(a[u]*a[v]));

      tau = LIBQP_MIN(LIBQP_MAX(tau,tau_lb),tau_ub);

      x[u] += tau/a[u];
      x[v] -= tau/a[v];

      /* update Nabla */
      for(i = 0; i < n; i++ ) 
         Nabla[i] += col_u[i]*tau/a[u] - col_v[i]*tau/a[v];

    }

    if( state.nIter >= MaxIter )
      state.exitflag = 0;

    if( print_state != NULL) 
    {
      state.QP = 0;
      for(i = 0; i < n; i++ ) 
        state.QP += 0.5*(x[i]*Nabla[i]+x[i]*f[i]); 

      print_state( state );
    }

  }  

  /* compute primal objective value */
  state.QP = 0;
  for(i = 0; i < n; i++ ) 
    state.QP += 0.5*(x[i]*Nabla[i]+x[i]*f[i]); 

cleanup:  

  LIBQP_FREE(Nabla);

  return( state ); 
}



