/*=================================================================================
 * SVMOCAS_LBP Train linear SVM classifier for images represented by LBP features. 
 * 
 * Synopsis:
 *  [W,W0,stat]= svmocas_lbp(Images,imSize,Wins,winSize,height_of_pyramid,X0,y,C,Method,TolRel,TolAbs,QPBound,BufSize,MaxTime,verb) 
 *
 * Input:  
 *   Images [(im_H*im_W) x nImages (uint8)]
 *   imSize [2 x 1 (uint32)] imSize = [im_H im_W]
 *   Wins [4 x nExamples (uint32)]  [image_idx; top_left_col; top_left_row; mirror]
 *   winSize [2 x 1 (uint32)] [win_H win_W]
 *   height_of_pyramid [1 x 1 (double)]
 *   X0 [1 x 1 (double)]
 *   y [nExamples x 1 (double)] +1/-1
 *   C [1 x 1 (double)] OR [nExamples x 1 (double)]
 *   Method [1x1 (double)] 0 (BMRM) or 1 (OCAS)
 *   TolRel [1x1 (double)]
 *   TolAbs [1x1 (double)]
 *   QPBound [1x1 (double)]
 *   BufSize [1x1 (double)]
 *   MaxTime [1x1 (double)]
 *   verb [1x1 (bouble)]
 * Output:
 *   W [nDim x 1] Parameter vector
 *   W0 [1x1] Bias term
 *   stat [struct] 
 *
 * Copyright (C) 2008, 2009, 2010 Vojtech Franc, xfrancv@cmp.felk.cvut.cz
 *                                Soeren Sonnenburg, soeren.sonnenburg@first.fraunhofer.de
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
#include "ocas_lbp_helper.h"
#include "liblbp.h"

#define DEFAULT_METHOD 1
#define DEFAULT_TOLREL 0.01
#define DEFAULT_TOLABS 0.0
#define DEFAULT_QPVALUE 0.0
#define DEFAULT_BUFSIZE 500
#define DEFAULT_MAXTIME mxGetInf()
#define DEFAULT_VERB 1

/*======================================================================
  Main code plus interface to Matlab.
========================================================================*/

void mexFunction( int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[] )
{
  double C, TolRel, TolAbs, QPBound, trn_err, MaxTime;
  double *vec_C;   
  uint32_t num_of_Cs;
  uint32_t i, j, BufSize;
  uint16_t Method;
  int verb;
  ocas_return_value_T ocas;
  double *tmp;

  /* timing variables */
  double init_time;
  double total_time;

  total_time = get_time();
  init_time = total_time;

  if(nrhs < 8 || nrhs > 15)
     mexErrMsgTxt("Improper number of input arguments.\n\n"
                  "SVMOCAS_LBP train linear SVM classifier for images prepresented by LBP features. \n\n"
                  "Synopsis:\n"
                  "  [W,W0,stat]= svmocas_lbp(Images,imSize,Wins,winSize,height_of_pyramid,\n"
                  "              X0,y,C,Method,TolRel,TolAbs,QPBound,BufSize,MaxTime,verb) \n\n"
                  "Input:  \n"
                  "  Images [(im_H*im_W) x nImages (uint8)]\n"
                  "  imSize [2 x 1 (uint32)] imSize = [im_H im_W]\n"
                  "  Wins [4 x nExamples (uint32)] [img_idx; topleft_col; topleft_row; mirror]\n"
                  "  winSize [2 x 1 (uint32)] [win_H win_W]\n"
                  "  height_of_pyramid [1 x 1 (double)]\n"
                  "  X0 [1 x 1 (double)]\n"
                  "  y [nExamples x 1 (double)] +1 or -1\n"
                  "  C [1 x 1 (double)] OR [nExamples x 1 (double)]\n"
                  "  Method [1x1 (double)] 0 for BMRM; 1 for OCAS \n"
                  "  TolRel [1x1 (double)]\n"
                  "  TolAbs [1x1 (double)]\n"
                  "  QPBound [1x1 (double)]\n"
                  "  BufSize [1x1 (double)]\n"
                  "  MaxTime [1x1 (double)]\n"
                  "  verb [1x1 (bouble)]\n\n"
                  "Output:\n"
                  "  W [nDim x 1] Parameter vector\n"
                  "  W0 [1x1] Bias term\n"
                  "  stat [struct] \n");

  /*                             0      1     2     3         4     5  6 7   8      9      10      11     12     13       14     15*/
  /*  [W,W0,stat]= svmocas_lbp(Images,imSize,Wins,winSize,nPyramids,X0,y,C,Method,TolRel,TolAbs,QPBound,BufSize,nExamples,MaxTime,verb) */

  if(nrhs >= 16)
    verb = (int)mxGetScalar(prhs[15]);
  else
    verb = DEFAULT_VERB;

  Images = (uint8_t*)mxGetPr(prhs[0]);
  nImages = mxGetN(prhs[0]);

  tmp = (double*)mxGetPr(prhs[1]);
  im_H = (uint32_t)tmp[0];
  im_W = (uint32_t)tmp[1];

  mexPrintf("im_h=%d  im_W=%d \n", im_H, im_W);
  if(mxGetM(prhs[0]) != im_H*im_W)
    mexErrMsgTxt("Dimension of Images does not match to im_H*im_W.");

  Wins = (uint32_t*)mxGetPr(prhs[2]);

  tmp = (double*)mxGetPr(prhs[3]);
  win_H = (uint32_t)tmp[0];
  win_W = (uint32_t)tmp[1];

  nPyramids = (uint32_t)mxGetScalar(prhs[4]);
/*  nDim = lbppyr_get_dim(win_H,win_W,nPyramids);*/
  nDim = liblbp_pyr_get_dim(win_H,win_W,nPyramids);

  croped_window = (uint32_t*)mxCalloc(win_H*win_W,sizeof(uint32_t));
  if(croped_window == NULL) 
    mexErrMsgTxt("Not enough memory for croped_window.");
  
  X0 = mxGetScalar(prhs[5]);
  data_y = (double*)mxGetPr(prhs[6]);

  nData = LIBOCAS_MAX(mxGetM(prhs[6]),mxGetN(prhs[6]));
  if(nData != mxGetN(prhs[2]))
    mexErrMsgTxt("Dimension missmatch betwenn Wins and y.");

  if(verb)
  {
    mexPrintf("Input data:\n"
              "   # of images     : %d\n"
              "   image height    : %d\n"
              "   image width     : %d\n",
              nImages, im_H, im_W);

    mexPrintf("Feature represenation:\n"
              "   base window height        : %d\n"
              "   base window width         : %d\n"
              "   nPyramids                 : %d\n"
              "   # of virtual examples     : %d\n"
              "   # of features per example : %d\n",
              win_H, win_W, nPyramids, nData, nDim);
  }

  num_of_Cs = LIBOCAS_MAX(mxGetN(prhs[7]),mxGetM(prhs[7]));
  if(num_of_Cs == 1)
  {
    C = (double)mxGetScalar(prhs[7]);
  }
  else
  {
    vec_C = (double*)mxGetPr(prhs[7]);
  }

  if(nrhs >= 9)
    Method = (uint32_t)mxGetScalar(prhs[8]);
  else
    Method = DEFAULT_METHOD;

  if(nrhs >= 10)
    TolRel = (double)mxGetScalar(prhs[9]);
  else
    TolRel = DEFAULT_TOLREL;
  
  if(nrhs >= 11)    
    TolAbs = (double)mxGetScalar(prhs[10]);
  else
    TolAbs = DEFAULT_TOLABS;

  if(nrhs >= 12)
    QPBound = (double)mxGetScalar(prhs[11]);
  else
    QPBound = DEFAULT_QPVALUE;

  if(nrhs >= 13)
    BufSize = (uint32_t)mxGetScalar(prhs[12]);
  else
    BufSize = DEFAULT_BUFSIZE;

  if(num_of_Cs > 1 && num_of_Cs < nData)
    mexErrMsgTxt("Length of the vector C less than the number of examples.");

  if(nrhs >= 14 && !mxIsInf(mxGetScalar(prhs[13])))
  {
    if((uint32_t)mxGetScalar(prhs[13]) < 0 || (uint32_t)mxGetScalar(prhs[13]) > nData)
      mexErrMsgTxt("Improper number of examples; must be > 0 and < max number of virtual example.\n");

    nData = (uint32_t)mxGetScalar(prhs[13]);
    mexPrintf("   # of examples set to : %d\n",nData);
  }

  if(nrhs >= 15)
    MaxTime = (double)mxGetScalar(prhs[14]);
  else
    MaxTime = DEFAULT_MAXTIME;

  /*----------------------------------------------------------------
    Print setting
  -------------------------------------------------------------------*/
  if(verb)
  {
    mexPrintf("SVM setting:\n");

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
/*  new_a = (double*)mxCalloc(nDim,sizeof(double));*/
  new_a = mxCalloc(nDim,sizeof(new_a[0]));
  if(new_a == NULL) 
    mexErrMsgTxt("Not enough memory for auxciliary cutting plane buffer new_a.");  

  if(num_of_Cs > 1)
  {
    for(i=0; i < nData; i++) 
      data_y[i] = data_y[i]*vec_C[i];
  }

  /* !!!!!!!!!!!!
  ptr = mxGetPr(data_X);
  for(i=0; i < nData; i++) {
    for(j=0; j < nDim; j++ ) {
      ptr[LIBOCAS_INDEX(j,i,nDim)] = ptr[LIBOCAS_INDEX(j,i,nDim)]*data_y[i];
    }
  }
  */

  /* init cutting plane buffer */
/*  full_A = mxCalloc(BufSize*nDim,sizeof(double));*/
  full_A = mxCalloc(BufSize*nDim,sizeof(full_A[0]));
  if( full_A == NULL )
    mexErrMsgTxt("Not enough memory for cutting plane buffer full_A."); 


  if(verb)
  {
    mexPrintf("Memory occupancy:\n"
             "   raw images         : %.2f MB\n" 
             "   CP buffer          : %.2f MB\n"
             "   parameter vector W : %.2f MB\n",
              (double)nImages*im_H*im_W/(1024*1024),
              (double)sizeof(full_A[0])*BufSize*nDim/(1024*1024), 
              (double)sizeof(W[0])*nDim/(1024*1024));
  }

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

  if(num_of_Cs > 1)
  {
    for(i=0; i < nData; i++) 
      data_y[i] = data_y[i]/vec_C[i];
  }

  plhs[1] = mxCreateDoubleScalar( W0 );
  
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

