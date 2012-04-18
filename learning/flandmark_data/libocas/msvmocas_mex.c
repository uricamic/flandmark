/*=================================================================================
 * svmocas_mex.c: Matlab MEX interface for OCAS solver training the linear SVM classifiers.
 * 
 * Synopsis:
 *  [W,stat] = msvmocas(X,y,C,Method,TolRel,TolAbs,QPBound,BufSize,nData,MaxTime,verb)
 *  [W,stat] = msvmocas(svmlight_data_file,C,Method,TolRel,TolAbs,QPBound,BufSize,nData,MaxTime,verb)
 *
 * Synopsis:
 *  [W,stat] = msvmocas(X,y,C,Method,TolRel,TolAbs,QPBound,BufSize,nData,MaxTime)
 * Copyright (C) 2008 Vojtech Franc, xfrancv@cmp.felk.cvut.cz
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public 
 * License as published by the Free Software Foundation; 
 *======================================================================================*/ 

#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <mex.h>

/*#define LIBOCAS_MATLAB*/

#include "libocas.h"
#include "ocas_helper.h"

#define DEFAULT_METHOD 1
#define DEFAULT_TOLREL 0.01
#define DEFAULT_TOLABS 0.0
#define DEFAULT_QPVALUE 0.0
#define DEFAULT_BUFSIZE 2000
#define DEFAULT_MAXTIME mxGetInf()
#define DEFAULT_VERB 1



/*======================================================================
  Main code plus interface to Matlab.
========================================================================*/

void mexFunction( int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[] )
{
  double C, TolRel, TolAbs, MaxTime, trn_err, QPBound;
  double *ptr;
  uint32_t i, j, BufSize;
  uint16_t Method;
  int verb;
  ocas_return_value_T ocas;

  /* timing variables */
  double init_time;
  double total_time;

  total_time = get_time();
  init_time = total_time;

  if(nrhs < 1)
    mexErrMsgTxt("Improper number of input arguments.");

  /* get input arguments */ 
  if(mxIsChar(prhs[0]) == false) 
  {

    if(nrhs < 3 || nrhs > 11)
      mexErrMsgTxt("Improper number of input arguments.");

    /*  [W,stat] = msvmocas(X,y,C,Method,TolRel,TolAbs,QPBound,BufSize,nData,MaxTime)*/

    data_X = (mxArray*)prhs[0];
    if (!(mxIsDouble(data_X)))
      mexErrMsgTxt("Input argument X must be of type double.");

    if (mxGetNumberOfDimensions(data_X) != 2)
      mexErrMsgTxt("Input argument X must be two dimensional.");

    data_y = (double*)mxGetPr(prhs[1]);

    if(LIBOCAS_MAX(mxGetM(prhs[1]),mxGetN(prhs[1])) != mxGetN(prhs[0]))
      mexErrMsgTxt("Length of vector y must equal to the number of columns of matrix X.");

    C = (double)mxGetScalar(prhs[2]);

    if(nrhs >= 4)
      Method = (uint32_t)mxGetScalar(prhs[3]);
    else
      Method = DEFAULT_METHOD;

    if(nrhs >= 5)
      TolRel = (double)mxGetScalar(prhs[4]);
    else
      TolRel = DEFAULT_TOLREL;

    if(nrhs >= 6)    
      TolAbs = (double)mxGetScalar(prhs[5]);
    else
      TolAbs = DEFAULT_TOLABS;

    if(nrhs >= 7)
      QPBound = (double)mxGetScalar(prhs[6]);
    else
      QPBound = DEFAULT_QPVALUE;
    
    if(nrhs >= 8)
      BufSize = (uint32_t)mxGetScalar(prhs[7]);
    else
      BufSize = DEFAULT_BUFSIZE;

    if(nrhs >= 9 && mxIsInf(mxGetScalar(prhs[8])) == false)
      nData = (uint32_t)mxGetScalar(prhs[8]);
    else
      nData = mxGetN(data_X);
      
    if(nData < 1 || nData > mxGetN(prhs[0])) 
      mexErrMsgTxt("Improper value of argument nData.");

    if(nrhs >= 10)
      MaxTime = (double)mxGetScalar(prhs[9]);
    else
      MaxTime = DEFAULT_MAXTIME;

    if(nrhs >= 11)
      verb = (int)mxGetScalar(prhs[10]);
    else
      verb = DEFAULT_VERB;


  }
  else
  {
    /*  [W,stat] = msvmocas(svmlight_data_file,C,Method,TolRel,TolAbs,QPBound,BufSize,nData,MaxTime)*/
    char *fname;
    int fname_len;

    if(nrhs < 2 || nrhs > 10)
      mexErrMsgTxt("Improper number of input arguments.");

    if(!mxIsChar(prhs[0]))
      mexErrMsgTxt("First input argument must be of type string.");

    fname_len = mxGetNumberOfElements(prhs[0]) + 1;   
    fname = mxCalloc(fname_len, sizeof(char));    

    if (mxGetString(prhs[0], fname, fname_len) != 0)     
      mexErrMsgTxt("Could not convert first input argument to string.");

    if(nrhs >= 10)
      verb = (int)mxGetScalar(prhs[9]);
    else
      verb = DEFAULT_VERB;

    /* load data */
    if( load_svmlight_file(fname,verb) == -1 || data_X == NULL || data_y == NULL)
      mexErrMsgTxt("Cannot load input file.");

    C = (double)mxGetScalar(prhs[1]);

    if(nrhs >= 3)
      Method = (uint32_t)mxGetScalar(prhs[2]);
    else
      Method = DEFAULT_METHOD;

    if(nrhs >= 4)
      TolRel = (double)mxGetScalar(prhs[3]);
    else
      TolRel = DEFAULT_TOLREL;

    if(nrhs >= 5)    
      TolAbs = (double)mxGetScalar(prhs[4]);
    else
      TolAbs = DEFAULT_TOLABS;

    if(nrhs >= 6)
      QPBound = (double)mxGetScalar(prhs[5]);
    else
      QPBound = DEFAULT_QPVALUE;
    
    if(nrhs >= 7)
      BufSize = (uint32_t)mxGetScalar(prhs[6]);
    else
      BufSize = DEFAULT_BUFSIZE;

    if(nrhs >= 8 && mxIsInf(mxGetScalar(prhs[7])) == false)
      nData = (uint32_t)mxGetScalar(prhs[7]);
    else
      nData = mxGetN(data_X);

    if(nData < 1 || nData > mxGetN(data_X)) 
      mexErrMsgTxt("Improper value of argument nData.");

    if(nrhs >= 9)
      MaxTime = (double)mxGetScalar(prhs[8]);
    else
      MaxTime = DEFAULT_MAXTIME;


  }

/*  nDim = mxGetM(prhs[0]);*/
  nDim = mxGetM(data_X);
  for(i=0, nY = 0; i < nData; i++)
    nY = LIBOCAS_MAX(nY, (uint32_t)data_y[i]);


  /*----------------------------------------------------------------
    Print setting
  -------------------------------------------------------------------*/
  if(verb)
  {
    mexPrintf("Input data statistics:\n"
              "   # of examples  : %d\n"
              "   # of classes   : %d\n"
              "   dimensionality : %d\n",
              nData, nY, nDim);
    
    if( mxIsSparse(data_X)== true ) 
      mexPrintf("   density        : %.2f%%\n",
                100.0*(double)mxGetNzmax(data_X)/((double)nDim*(double)(mxGetN(data_X))));
    else
      mexPrintf("    density       : 100%% (full)\n");

    mexPrintf("Setting:\n"
         "   C              : %f\n"
         "   # of examples  : %d\n"
         "   solver         : %d\n"
         "   cache size     : %d\n"
         "   TolAbs         : %f\n"
         "   TolRel         : %f\n"
         "   QPValue        : %f\n"
         "   MaxTime        : %f [s]\n",
         C, nData, Method,BufSize,TolAbs,TolRel, QPBound, MaxTime);
  }
  
  /* learned weight vector */
  plhs[0] = (mxArray*)mxCreateDoubleMatrix(nDim,nY,mxREAL);
  W = (double*)mxGetPr(plhs[0]);
  if(W == NULL) mexErrMsgTxt("Not enough memory for vector W.");

  oldW = (double*)mxCalloc(nY*nDim,sizeof(double));
  if(oldW == NULL) mexErrMsgTxt("Not enough memory for vector oldW.");

  /* allocate buffer for computing cutting plane */
  new_a = (double*)mxCalloc(nY*nDim,sizeof(double));
  if(new_a == NULL) 
    mexErrMsgTxt("Not enough memory for auxciliary cutting plane buffer new_a.");  


  if( mxIsSparse(data_X)== true ) 
  {
    /* init cutting plane buffer */
    sparse_A.nz_dims = mxCalloc(BufSize,sizeof(uint32_t));
    sparse_A.index = mxCalloc(BufSize,sizeof(sparse_A.index[0]));
    sparse_A.value = mxCalloc(BufSize,sizeof(sparse_A.value[0]));
    if(sparse_A.nz_dims == NULL || sparse_A.index == NULL || sparse_A.value == NULL) 
      mexErrMsgTxt("Not enough memory for cutting plane buffer sparse_A.");  

    if(verb)
      mexPrintf("Starting optimization:\n");

    init_time=get_time()-init_time;

    if(verb)
      ocas = msvm_ocas_solver( C, data_y, nY, nData, TolRel, TolAbs, QPBound, MaxTime,BufSize, Method,
                             &msvm_sparse_compute_W, &msvm_full_update_W, &msvm_sparse_add_new_cut,
                             &msvm_sparse_compute_output, &qsort_data, &ocas_print, 0);
    else
      ocas = msvm_ocas_solver( C, data_y, nY, nData, TolRel, TolAbs, QPBound, MaxTime,BufSize, Method,
                             &msvm_sparse_compute_W, &msvm_full_update_W, &msvm_sparse_add_new_cut,
                             &msvm_sparse_compute_output, &qsort_data, &ocas_print_null, 0);
  }
  else
  {
    /* init cutting plane buffer */
    full_A = mxCalloc(BufSize*nDim*nY,sizeof(double));
    if( full_A == NULL )
      mexErrMsgTxt("Not enough memory for cutting plane buffer full_A.");    

    if(verb)
      mexPrintf("Starting optimization:\n");

    init_time=get_time()-init_time;

    if(verb)
      ocas = msvm_ocas_solver( C, data_y, nY, nData, TolRel, TolAbs, QPBound, MaxTime,BufSize, Method,
                             &msvm_full_compute_W, &msvm_full_update_W, &msvm_full_add_new_cut,
                             &msvm_full_compute_output, &qsort_data, &ocas_print, 0); 
    else
      ocas = msvm_ocas_solver( C, data_y, nY, nData, TolRel, TolAbs, QPBound, MaxTime,BufSize, Method,
                             &msvm_full_compute_W, &msvm_full_update_W, &msvm_full_add_new_cut,
                             &msvm_full_compute_output, &qsort_data, &ocas_print_null, 0); 

  }
 
  if(verb)
  {
    mexPrintf("Stopping condition: ");
    switch( ocas.exitflag )
    {
       case 1: mexPrintf("1-Q_D/Q_P <= TolRel(=%f) satisfied.\n", TolRel); break;
       case 2: mexPrintf("Q_P-Q_D <= TolAbs(=%f) satisfied.\n", TolAbs); break;
       case 3: mexPrintf("Q_P <= QPBound(=%f) satisfied.\n", QPBound); break;
       case 4: mexPrintf("Optimization time (=%f) >= MaxTime(=%f).\n", ocas.ocas_time, MaxTime); break;
       case -1: mexPrintf("Has not converged!\n" ); break;
       case -2: mexPrintf("Not enough memory for the solver.\n" ); break;
    }
  }

  total_time=get_time()-total_time;
  if(verb)
  {
    mexPrintf("Timing statistics:\n"
			"   init_time      : %f[s]\n"
			"   qp_solver_time : %f[s]\n"
			"   sort_time      : %f[s]\n"
			"   output_time    : %f[s]\n"
			"   add_time       : %f[s]\n"
			"   w_time         : %f[s]\n"
			"   print_time     : %f[s]\n"
			"   ocas_time      : %f[s]\n"
			"   total_time     : %f[s]\n",
			init_time, ocas.qp_solver_time, ocas.sort_time, ocas.output_time, 
            ocas.add_time, ocas.w_time, ocas.print_time, ocas.ocas_time, total_time);

    mexPrintf("Training error: %.4f%%\n", 100*(double)ocas.trn_err/(double)nData);
  }

  
  /* return ocas optimizer statistics */
  /* typedef struct {
     uint32_t nIter;    
     uint32_t nCutPlanes;
     double trn_err;      
     double Q_P;          
     double Q_D;
     double output_time;
     double sort_time;
     double solver_time;
     int8_t exitflag;       
     } ocas_return_value_T; */

  const char *field_names[] = {"nTrnErrors","Q_P","Q_D","nIter","nCutPlanes","exitflag",
                               "init_time","output_time","sort_time","qp_solver_time","add_time",
                               "w_time","ocas_time","total_time"}; 
  mwSize dims[2] = {1,1};  

  plhs[1] = mxCreateStructArray(2, dims, (sizeof(field_names)/sizeof(*field_names)), field_names);
  
  mxSetField(plhs[1],0,"nIter",mxCreateDoubleScalar((double)ocas.nIter));
  mxSetField(plhs[1],0,"nCutPlanes",mxCreateDoubleScalar((double)ocas.nCutPlanes));
  mxSetField(plhs[1],0,"nTrnErrors",mxCreateDoubleScalar(ocas.trn_err)); 
  mxSetField(plhs[1],0,"Q_P",mxCreateDoubleScalar(ocas.Q_P)); 
  mxSetField(plhs[1],0,"Q_D",mxCreateDoubleScalar(ocas.Q_D)); 
  mxSetField(plhs[1],0,"init_time",mxCreateDoubleScalar(init_time)); 
  mxSetField(plhs[1],0,"output_time",mxCreateDoubleScalar(ocas.output_time)); 
  mxSetField(plhs[1],0,"sort_time",mxCreateDoubleScalar(ocas.sort_time)); 
  mxSetField(plhs[1],0,"qp_solver_time",mxCreateDoubleScalar(ocas.qp_solver_time)); 
  mxSetField(plhs[1],0,"add_time",mxCreateDoubleScalar(ocas.add_time)); 
  mxSetField(plhs[1],0,"w_time",mxCreateDoubleScalar(ocas.w_time)); 
  mxSetField(plhs[1],0,"ocas_time",mxCreateDoubleScalar(ocas.ocas_time)); 
  mxSetField(plhs[1],0,"total_time",mxCreateDoubleScalar(total_time)); 
  mxSetField(plhs[1],0,"exitflag",mxCreateDoubleScalar((double)ocas.exitflag)); 

  return;
}

