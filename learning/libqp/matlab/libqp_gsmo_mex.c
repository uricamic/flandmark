/*--------------------------------------------------------------------------
 gsmo_mex.c: Matlab MEX interface for the Generalized SMO solver.

 Synopsis:
  [x,QP,exitflag,nIter] = ...
     libqp_gsmo_mex(H,f,a,b,LB,UB,x0,MaxIter,TolKKT,verb)

 Compile:
  mex gsmo_mex.c gsmolib.c

 Description:
   See "help libqp_gsmo"
                                                                    
 Copyright (C) 2006-2008 Vojtech Franc, xfrancv@cmp.felk.cvut.cz
 Center for Machine Perception, CTU FEL Prague

-------------------------------------------------------------------- */

#include "mex.h"
#include "string.h"

#define LIBQP_MATLAB
#include "../lib/libqp.h"


#define INDEX(ROW,COL,NUM_ROWS) ((COL*NUM_ROWS)+ROW)
#define MIN(A,B) ((A < B) ? A : B)
#define MAX(A,B) ((A > B) ? A : B)

/* ------------------------------------------------------------*/
/* Declaration of global variables                             */
/* ------------------------------------------------------------*/
double *mat_H;
uint32_t nVar;

/* ------------------------------------------------------------
  Returns pointer at the a-th column of the matrix H.
------------------------------------------------------------ */
const double *get_col( uint32_t i )
{
  return( &mat_H[ nVar*i ] );
}


/* -------------------------------------------------------------------
 Main MEX function - interface to Matlab.

  [vec_x,exitflag,t,access,Nabla] = 
         gsmo_mex(H,f,a,b,LB,UB,x0,Nabla0,tmax,tolKKT,verb);

-------------------------------------------------------------------- */
void mexFunction( int nlhs, mxArray *plhs[],int nrhs, const mxArray*prhs[] )
{
  int verb;          
  uint32_t i;         
  uint32_t MaxIter;   
  double TolKKT; 
  double *vec_x;         /* output arg -- solution*/ 
  double *vec_x0;         
  double *diag_H;    /* diagonal of matrix H */
  double *f;         /* vector f */
  double *a;
  double b;
  double *LB;
  double *UB;
  double fval;
  
  /*------------------------------------------------------------------- */
  /* Take input arguments                                               */
  /*------------------------------------------------------------------- */

  if( nrhs != 10) mexErrMsgTxt("Incorrect number of input arguments.");

  mat_H = mxGetPr(prhs[0]);
  nVar = mxGetM(prhs[0]);
  f = mxGetPr(prhs[1]);   
  a = mxGetPr(prhs[2]);
  b = mxGetScalar(prhs[3]);
  LB = mxGetPr(prhs[4]);
  UB = mxGetPr(prhs[5]);
  vec_x0 = mxGetPr(prhs[6]);
  MaxIter = mxIsInf( mxGetScalar(prhs[7])) ? INT_MAX : (long)mxGetScalar(prhs[7]);
  TolKKT = mxGetScalar(prhs[8]); 
  verb = (int)(mxGetScalar(prhs[9])); 
    
  if( verb > 0 ) {
    mexPrintf("Settings of QP solver\n");
    mexPrintf("MaxIter : %d\n", MaxIter );
    mexPrintf("TolKKT  : %f\n", TolKKT );
    mexPrintf("nVar    : %d\n", nVar );
    mexPrintf("verb    : %d\n", verb );
  }

  plhs[0] = mxCreateDoubleMatrix(nVar,1,mxREAL);
  vec_x = mxGetPr(plhs[0]);
  memcpy( vec_x, vec_x0, sizeof(double)*nVar );

  diag_H = mxCalloc(nVar, sizeof(double));
  if( diag_H == NULL ) mexErrMsgTxt("Not enough memory.");
  for(i = 0; i < nVar; i++ ) 
    diag_H[i] = mat_H[nVar*i+i];
  

  /*------------------------------------------------------------------- */
  /* Call QP solver                                                     */
  /*------------------------------------------------------------------- */
  libqp_state_T state;
  state = libqp_gsmo_solver( &get_col, diag_H, f, a, b, LB, UB, vec_x, 
                        nVar, MaxIter, TolKKT, NULL );

  /*------------------------------------------------------------------- */
  /* Generate outputs                                                   */
  /*------------------------------------------------------------------- */

  plhs[1] = mxCreateDoubleScalar(state.QP);
  plhs[2] = mxCreateDoubleScalar((double)state.exitflag);
  plhs[3] = mxCreateDoubleScalar((double)state.nIter);


  /*------------------------------------------------------------------- */
  /* Clean up                                                           */
  /*------------------------------------------------------------------- */
  mxFree( diag_H );
}

/* EOF */
