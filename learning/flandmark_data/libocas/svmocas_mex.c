/*=================================================================================
 * svmocas_mex.c: Matlab MEX interface for OCAS solver training the linear SVM classifiers.
 *
 * Synopsis:
 *  [W,W0,stat] = svmocas(X,X0,y,C,Method,TolRel,TolAbs,QPBound,BufSize,nData,MaxTime,verb)
 *  [W,W0,stat] = svmocas(data_file,X0,C,Method,TolRel,TolAbs,QPBound,BufSize,nData,MaxTime,verb)
 *
 * See svmocas.m for more help.
 *
 * Copyright (C) 2008, 2009 Vojtech Franc, xfrancv@cmp.felk.cvut.cz
 *                          Soeren Sonnenburg, soeren.sonnenburg@first.fraunhofer.de
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public 
 * License as published by the Free Software Foundation; 
 *======================================================================================*/ 

#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <mex.h>

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
  double C, TolRel, TolAbs, QPBound, trn_err, MaxTime;
  double *vec_C;   
  double *ptr;
  uint32_t num_of_Cs;
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
    /* [W,W0,stat] = svmocas_mex(X,X0,y,C,Method,TolRel,TolAbs,QPBound,BufSize,nData,MaxTime); */

    if(nrhs < 4 || nrhs > 12)
      mexErrMsgTxt("Improper number of input arguments.");

    if(nrhs >= 12)
      verb = (int)mxGetScalar(prhs[11]);
    else
      verb = DEFAULT_VERB;


    data_X = (mxArray*)prhs[0];
    if (!(mxIsDouble(data_X)))
        mexErrMsgTxt("Input argument X must be of type double.");

    if (mxGetNumberOfDimensions(data_X) != 2)
        mexErrMsgTxt("Input argument X must be two dimensional");

    X0 = mxGetScalar(prhs[1]);
    data_y = (double*)mxGetPr(prhs[2]);

    if(LIBOCAS_MAX(mxGetM(prhs[2]),mxGetN(prhs[2])) != mxGetN(prhs[0]))
      mexErrMsgTxt("Length of vector y must equl to the number of columns of matrix X.");

    nDim = mxGetM(prhs[0]);

    if(verb)
    {
      mexPrintf("Input data statistics:\n"
                "   # of examples  : %d\n"
                "   dimensionality : %d\n",
                mxGetN(data_X), nDim);
    
      if( mxIsSparse(data_X)== true ) 
        mexPrintf("   density        : %.2f%%\n",
                  100.0*(double)mxGetNzmax(data_X)/((double)nDim*(double)(mxGetN(data_X))));
      else
        mexPrintf("   density        : 100%% (full)\n");
    }


    num_of_Cs = LIBOCAS_MAX(mxGetN(prhs[3]),mxGetM(prhs[3]));
    if(num_of_Cs == 1)
    {
       C = (double)mxGetScalar(prhs[3]);
    }
    else
    {
       vec_C = (double*)mxGetPr(prhs[3]);
    }

    if(nrhs >= 5)
      Method = (uint32_t)mxGetScalar(prhs[4]);
    else
      Method = DEFAULT_METHOD;

    if(nrhs >= 6)
      TolRel = (double)mxGetScalar(prhs[5]);
    else
      TolRel = DEFAULT_TOLREL;

    if(nrhs >= 7)    
      TolAbs = (double)mxGetScalar(prhs[6]);
    else
      TolAbs = DEFAULT_TOLABS;

    if(nrhs >= 8)
      QPBound = (double)mxGetScalar(prhs[7]);
    else
      QPBound = DEFAULT_QPVALUE;

    if(nrhs >= 9)
      BufSize = (uint32_t)mxGetScalar(prhs[8]);
    else
      BufSize = DEFAULT_BUFSIZE;

    if(nrhs >= 10 && mxIsInf(mxGetScalar(prhs[9])) == false)
      nData = (uint32_t)mxGetScalar(prhs[9]);
    else
      nData = mxGetN(data_X);
      
    if(nData < 1 || nData > mxGetN(prhs[0])) 
      mexErrMsgTxt("Improper value of argument nData.");

    if(num_of_Cs > 1 && num_of_Cs < nData)
      mexErrMsgTxt("Length of the vector C less than the number of examples.");

    if(nrhs >= 11)
      MaxTime = (double)mxGetScalar(prhs[10]);
    else
      MaxTime = DEFAULT_MAXTIME;
  } 
  else
  {
    /* [W,W0,stat] = svmocas_mex(svmlight_data_file,X0,C,Method,TolRel,TolAbs,QPBound,BufSize,nData,MaxTime); */
    char *fname;
    int fname_len;

    if(nrhs < 3 || nrhs > 11)
      mexErrMsgTxt("Improper number of input arguments.");

    if(nrhs >= 11)
      verb = (int)mxGetScalar(prhs[10]);
    else
      verb = DEFAULT_VERB;

    if(!mxIsChar(prhs[0]))
      mexErrMsgTxt("First input argument must be of type string.");

    fname_len = mxGetNumberOfElements(prhs[0]) + 1;   
    fname = mxCalloc(fname_len, sizeof(char));    

    if (mxGetString(prhs[0], fname, fname_len) != 0)     
      mexErrMsgTxt("Could not convert first input argument to string.");

    if( load_svmlight_file(fname,verb) == -1 || data_X == NULL || data_y == NULL)
      mexErrMsgTxt("Cannot load input file.");

    nDim = mxGetM(data_X);
    X0 = mxGetScalar(prhs[1]);

/*    C = (double)mxGetScalar(prhs[2]);*/
    num_of_Cs = LIBOCAS_MAX(mxGetN(prhs[2]),mxGetM(prhs[2]));
    if(num_of_Cs == 1)
    {
       C = (double)mxGetScalar(prhs[2]);
    }
    else
    {
       vec_C = (double*)mxGetPr(prhs[2]);
    }

    if(verb)
      mexPrintf("Input data statistics:\n"
                "   # of examples  : %d\n"
                "   dimensionality : %d\n"
                "   density        : %.2f%%\n",
                mxGetN(data_X), nDim, 100.0*(double)mxGetNzmax(data_X)/((double)nDim*(double)(mxGetN(data_X))));
    
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

    if(nData < 1 || nData > mxGetN(data_X)) 
      mexErrMsgTxt("Improper value of argument nData.");

    if(num_of_Cs > 1 && num_of_Cs < nData)
      mexErrMsgTxt("Length of the vector C less than the number of examples.");

    if(nrhs >= 10)
      MaxTime = (double)mxGetScalar(prhs[9]);
    else
      MaxTime = DEFAULT_MAXTIME;

  } 


  /*----------------------------------------------------------------
    Print setting
  -------------------------------------------------------------------*/
  if(verb)
  {
    mexPrintf("Setting:\n");

    if( num_of_Cs == 1)
      mexPrintf("   C              : %f\n", C);
    else
      mexPrintf("   C              : different for each example\n");

    mexPrintf("   bias           : %.0f\n"
              "   # of examples  : %d\n"
              "   solver         : %d\n"
              "   cache size     : %d\n"
              "   TolAbs         : %f\n"
              "   TolRel         : %f\n"
              "   QPValue        : %f\n"
              "   MaxTime        : %f [s]\n"
              "   verb           : %d\n",
              X0, nData, Method,BufSize,TolAbs,TolRel, QPBound, MaxTime, verb);
  }
  
  /* learned weight vector */
  plhs[0] = (mxArray*)mxCreateDoubleMatrix(nDim,1,mxREAL);
  W = (double*)mxGetPr(plhs[0]);
  if(W == NULL) mexErrMsgTxt("Not enough memory for vector W.");

  oldW = (double*)mxCalloc(nDim,sizeof(double));
  if(oldW == NULL) mexErrMsgTxt("Not enough memory for vector oldW.");

  W0 = 0;
  oldW0 = 0;

  A0 = mxCalloc(BufSize,sizeof(A0[0]));
  if(A0 == NULL) mexErrMsgTxt("Not enough memory for vector A0.");

  /* allocate buffer for computing cutting plane */
  new_a = (double*)mxCalloc(nDim,sizeof(double));
  if(new_a == NULL) 
    mexErrMsgTxt("Not enough memory for auxciliary cutting plane buffer new_a.");  


  if(num_of_Cs > 1)
  {
    for(i=0; i < nData; i++) 
      data_y[i] = data_y[i]*vec_C[i];
  }


  if( mxIsSparse(data_X)== true ) {

    /* for i=1:nData, X(:,i) = X(:,i)*y(i); end*/
    for(i=0; i < nData; i++) 
        mul_sparse_col(data_y[i], data_X, i);           
  

    /* init cutting plane buffer */
    sparse_A.nz_dims = mxCalloc(BufSize,sizeof(uint32_t));
    sparse_A.index = mxCalloc(BufSize,sizeof(sparse_A.index[0]));
    sparse_A.value = mxCalloc(BufSize,sizeof(sparse_A.value[0]));
    if(sparse_A.nz_dims == NULL || sparse_A.index == NULL || sparse_A.value == NULL) 
      mexErrMsgTxt("Not enough memory for cutting plane buffer sparse_A.");  

    init_time=get_time()-init_time;

    if(verb)
    {
      mexPrintf("Starting optimization:\n");

      if(num_of_Cs == 1)
      {
        ocas = svm_ocas_solver( C, nData, TolRel, TolAbs, QPBound, MaxTime,BufSize, Method, 
                              &sparse_compute_W, &sparse_update_W, &sparse_add_new_cut, 
                              &sparse_compute_output, &qsort_data, &ocas_print, 0);
      }  
      else
      {
        ocas = svm_ocas_solver_difC( vec_C, nData, TolRel, TolAbs, QPBound, MaxTime,BufSize, Method, 
                                     &sparse_compute_W, &sparse_update_W, &sparse_add_new_cut, 
                                     &sparse_compute_output, &qsort_data, &ocas_print, 0);
      }  

    }
    else
    {
      if(num_of_Cs == 1)
      {
        ocas = svm_ocas_solver( C, nData, TolRel, TolAbs, QPBound, MaxTime,BufSize, Method, 
                              &sparse_compute_W, &sparse_update_W, &sparse_add_new_cut, 
                              &sparse_compute_output, &qsort_data, &ocas_print_null, 0);
      }
      else
      {
        ocas = svm_ocas_solver_difC( vec_C, nData, TolRel, TolAbs, QPBound, MaxTime,BufSize, Method, 
                              &sparse_compute_W, &sparse_update_W, &sparse_add_new_cut, 
                              &sparse_compute_output, &qsort_data, &ocas_print_null, 0);
      }
    }

  }
  else
  {

    ptr = mxGetPr(data_X);
    for(i=0; i < nData; i++) {
      for(j=0; j < nDim; j++ ) {
        ptr[LIBOCAS_INDEX(j,i,nDim)] = ptr[LIBOCAS_INDEX(j,i,nDim)]*data_y[i];
      }
    }

    /* init cutting plane buffer */
    full_A = mxCalloc(BufSize*nDim,sizeof(double));
    if( full_A == NULL )
      mexErrMsgTxt("Not enough memory for cutting plane buffer full_A.");    

    init_time=get_time()-init_time;

    if(verb)
    {
      mexPrintf("Starting optimization:\n");
    
      if(num_of_Cs == 1)
        ocas = svm_ocas_solver( C, nData, TolRel, TolAbs, QPBound, MaxTime,BufSize, Method, 
                                &full_compute_W, &full_update_W, &full_add_new_cut, 
                                &full_compute_output, &qsort_data, &ocas_print, 0);
      else
        ocas = svm_ocas_solver_difC( vec_C, nData, TolRel, TolAbs, QPBound, MaxTime,BufSize, Method, 
                                &full_compute_W, &full_update_W, &full_add_new_cut, 
                                &full_compute_output, &qsort_data, &ocas_print, 0);

    }
    else
    {
      if(num_of_Cs == 1)
        ocas = svm_ocas_solver( C, nData, TolRel, TolAbs, QPBound, MaxTime,BufSize, Method, 
                                &full_compute_W, &full_update_W, &full_add_new_cut, 
                                &full_compute_output, &qsort_data, &ocas_print_null, 0);
      else
        ocas = svm_ocas_solver_difC( vec_C, nData, TolRel, TolAbs, QPBound, MaxTime,BufSize, Method, 
                                     &full_compute_W, &full_update_W, &full_add_new_cut, 
                                     &full_compute_output, &qsort_data, &ocas_print_null, 0);
    }
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

  /* multiply data ba labels as it was at the begining */
  if( mxIsSparse(data_X)== true ) {
    /* for i=1:nData, X(:,i) = X(:,i)*y(i); end*/
    for(i=0; i < nData; i++) 
    {
        mul_sparse_col(1/data_y[i], data_X, i);
    }
  }
  else
  {
    ptr = mxGetPr(data_X);
    for(i=0; i < nData; i++) {
      for(j=0; j < nDim; j++ ) {
        ptr[LIBOCAS_INDEX(j,i,nDim)] = ptr[LIBOCAS_INDEX(j,i,nDim)]/data_y[i];
      }
    }
  }

  if(num_of_Cs > 1)
  {
    for(i=0; i < nData; i++) 
      data_y[i] = data_y[i]/vec_C[i];
  }


  plhs[1] = mxCreateDoubleScalar( W0 );
  
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
                               "w_time","print_time","ocas_time","total_time"}; 
  mwSize dims[2] = {1,1};  

  plhs[2] = mxCreateStructArray(2, dims, (sizeof(field_names)/sizeof(*field_names)), field_names);
  
  mxSetField(plhs[2],0,"nIter",mxCreateDoubleScalar((double)ocas.nIter));
  mxSetField(plhs[2],0,"nCutPlanes",mxCreateDoubleScalar((double)ocas.nCutPlanes));
  mxSetField(plhs[2],0,"nTrnErrors",mxCreateDoubleScalar(ocas.trn_err)); 
  mxSetField(plhs[2],0,"Q_P",mxCreateDoubleScalar(ocas.Q_P)); 
  mxSetField(plhs[2],0,"Q_D",mxCreateDoubleScalar(ocas.Q_D)); 
  mxSetField(plhs[2],0,"init_time",mxCreateDoubleScalar(init_time)); 
  mxSetField(plhs[2],0,"output_time",mxCreateDoubleScalar(ocas.output_time)); 
  mxSetField(plhs[2],0,"sort_time",mxCreateDoubleScalar(ocas.sort_time)); 
  mxSetField(plhs[2],0,"qp_solver_time",mxCreateDoubleScalar(ocas.qp_solver_time)); 
  mxSetField(plhs[2],0,"add_time",mxCreateDoubleScalar(ocas.add_time)); 
  mxSetField(plhs[2],0,"w_time",mxCreateDoubleScalar(ocas.w_time)); 
  mxSetField(plhs[2],0,"print_time",mxCreateDoubleScalar(ocas.print_time)); 
  mxSetField(plhs[2],0,"ocas_time",mxCreateDoubleScalar(ocas.ocas_time)); 
  mxSetField(plhs[2],0,"total_time",mxCreateDoubleScalar(total_time)); 
  mxSetField(plhs[2],0,"exitflag",mxCreateDoubleScalar((double)ocas.exitflag)); 

  return;
}

