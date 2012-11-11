/* -----------------------------------------------------------------------
libqp_splx_mex.c: MEX-file interface for libqp_splx solver.
                                                                     
 Synopsis:                                                                     
  [x,QP,QD,exitflag,nIter] = libqp_splx_mex(H,diagH,f,b,I,S,x0,MaxIter,TolAbs,TolRel,QP_TH,verb)

 Compile:

  mex libqp_splx_mex.c libqp_splx.c

 Description:
   See "help libqp_splx".
                                                                    
 Copyright (C) 2006-2008 Vojtech Franc, xfrancv@cmp.felk.cvut.cz
 Center for Machine Perception, CTU FEL Prague

-------------------------------------------------------------------------*/

#include "mex.h"
#include "string.h"

#define LIBQP_MATLAB
#include "../lib/libqp.h"

/* -- Global variables --------------------------------------*/

mxArray  *mat_H;     /* pointer to the Hessian matrix [n x n] */
uint32_t nVar;       /* number of ariables */
double *col_buf0, *col_buf1;
int col_buf_idx = 0;


/* ------------------------------------------------------------
  Returns pointer at the i-th column of the Hessian matrix H.
------------------------------------------------------------ */
const double *get_col( uint32_t col )
{
  double *pr = mxGetPr(mat_H);

  /*  mexPrintf("col=%d\n", col);*/
  return( &pr[ nVar*col ] );
}

const double *get_col_sparse( uint32_t col )
{
  uint32_t nItems, ptr, i, row;
  mwIndex *Ir, *Jc;
  double *Pr, val;

  /*  mexPrintf("col=%d\n", col);*/

  Ir = mxGetIr(mat_H);
  Jc = mxGetJc(mat_H);
  Pr = mxGetPr(mat_H);

  nItems = Jc[col+1] - Jc[col];
  /*  mexPrintf("nItems=%d\n", nItems);*/
  ptr = Jc[col];
  
  /*  memset(column_of_H, 0, sizeof(double)*nVar);*/
  if(col_buf_idx == 0)
  {
    memset(col_buf0, 0, sizeof(double)*nVar);
    /*    for(i=0; i < nVar; i++) col_buf0[i] = 0;*/

    for(i=0; i < nItems; i++) 
    {
      val = Pr[ptr];
      row = Ir[ptr++];

      col_buf0[row] = val;
    }
    col_buf_idx = 1;
    return(col_buf0);
  }
  else
  {
    memset(col_buf1, 0, sizeof(double)*nVar);
    /*    for(i=0; i < nVar; i++) col_buf1[i] = 0;*/

    for(i=0; i < nItems; i++) 
    {
      val = Pr[ptr];
      row = Ir[ptr++];

      col_buf_idx = 0;
      col_buf1[row] = val;
    }
    return(col_buf1);
  }

  return(NULL);
}


/* -------------------------------------------------------------
  MEX main function.
------------------------------------------------------------- */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
  int verb;          
  uint32_t MaxIter;
  uint32_t *vec_I; 
  uint32_t nConstr; 
  uint8_t *vec_S;
  double *vec_x;         
  double TolRel;     
  double TolAbs;     
  double QP_TH;
  double *vec_x0;        
  double *vec_f;         
  double *vec_b;
  double *diag_H;    
  long i ;           
  libqp_state_T state;
  double *tmp;

  /*------------------------------------------------------------------- 
     Get input arguments
   ------------------------------------------------------------------- */
  if( nrhs != 12) mexErrMsgTxt("Incorrect number of input arguments.");
  
  mat_H = (mxArray*)prhs[0];
  nVar = mxGetM(prhs[0]);

  diag_H = mxGetPr(prhs[1]);

  vec_f = mxGetPr(prhs[2]);
  vec_b = (double*)mxGetPr(prhs[3]);
  vec_I = (uint32_t*)mxGetPr(prhs[4]);
  vec_S = (uint8_t*)mxGetPr(prhs[5]);

  nConstr = LIBQP_MAX(mxGetN(prhs[3]),mxGetM(prhs[3]));

  vec_x0 = mxGetPr(prhs[6]);
  MaxIter = mxIsInf( mxGetScalar(prhs[7])) ? 0xFFFFFFFF : (uint32_t)mxGetScalar(prhs[7]);
  TolAbs = mxGetScalar(prhs[8]);   
  TolRel = mxGetScalar(prhs[9]);   
  QP_TH = mxGetScalar(prhs[10]);   
  verb = (int)mxGetScalar(prhs[11]);  

  /* print input setting if required */  
  if( verb > 0 ) {
    mexPrintf("Settings of LIBQP_SSVM solver:\n");
    mexPrintf("MaxIter  : %u\n", MaxIter );
    mexPrintf("TolAbs   : %f\n", TolAbs );
    mexPrintf("TolRel   : %f\n", TolRel );
    mexPrintf("QP_TH    : %f\n", QP_TH );
    mexPrintf("nVar     : %u\n", nVar );
    mexPrintf("nConstr  : %u\n", nConstr );
  }     
  
  /*------------------------------------------------------------------- 
     Inicialization                                                     
   ------------------------------------------------------------------- */

  /* create solution vector x [nVar x 1] */
  plhs[0] = mxCreateDoubleMatrix(nVar,1,mxREAL);
  vec_x = mxGetPr(plhs[0]);
  memcpy( vec_x, vec_x0, sizeof(double)*nVar );
  
  /*------------------------------------------------------------------- 
   Call the QP solver.
   -------------------------------------------------------------------*/

  if( mxIsSparse(mat_H)== true )
  {
    col_buf0 = mxCalloc(nVar,sizeof(double));
    if( col_buf0 == NULL)
      mexErrMsgTxt("Cannot allocate memory for col_buf0.");      

    col_buf1 = mxCalloc(nVar,sizeof(double));
    if( col_buf1 == NULL)
      mexErrMsgTxt("Cannot allocate memory for col_buf1.");      

    /*    tmp = get_col_sparse(nVar-1);
    for(i=0; i<nVar; i++)
    {
    mexPrintf("%d: %f\n",i,tmp[i]);
      vec_x[i] = tmp[i];
    }
*/
    
    state = libqp_splx_solver(&get_col_sparse, diag_H, vec_f, vec_b, vec_I, vec_S, vec_x, nVar, 
           MaxIter, TolAbs, TolRel, QP_TH, NULL);

    mxFree(col_buf0);
    mxFree(col_buf1);
    
  }
  else
  {

    /*    tmp = get_col(nVar-1);
    for(i=0; i<nVar; i++)
    {
       mexPrintf("%d: %f\n",i,tmp[i]);
      vec_x[i] = tmp[i];
    }
    */
   
    state = libqp_splx_solver(&get_col, diag_H, vec_f, vec_b, vec_I, vec_S, vec_x, nVar, 
                              MaxIter, TolAbs, TolRel, QP_TH, NULL);
   
  }
  
  /*------------------------------------------------------------------- 
    Set output arguments                                                   
    [x,QP,QD,exitflag,nIter] = libqp_splx_mex(...)
  ------------------------------------------------------------------- */

  plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL);
  *(mxGetPr(plhs[1])) = (double)state.QP;

  plhs[2] = mxCreateDoubleMatrix(1,1,mxREAL);
  *(mxGetPr(plhs[2])) = (double)state.QD;

  plhs[3] = mxCreateDoubleMatrix(1,1,mxREAL);
  *(mxGetPr(plhs[3])) = (double)state.exitflag;

  plhs[4] = mxCreateDoubleMatrix(1,1,mxREAL);
  *(mxGetPr(plhs[4])) = (double)state.nIter;

  
  return; 
}
 
