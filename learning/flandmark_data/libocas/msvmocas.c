/*-----------------------------------------------------------------------
 * msvmocas.c: Standalone application implementing the OCAS folver for 
 *   training multi-class linear SVM classifiers.
 *   
 * Copyright (C) 2008,2009 Vojtech Franc, xfrancv@cmp.felk.cvut.cz
 *                    Soeren Sonnenburg, soeren.sonnenburg@first.fraunhofer.de
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public 
 * License as published by the Free Software Foundation; 
 *-------------------------------------------------------------------- */


#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <stdlib.h>

#include "lib_svmlight_format.h"
#include "libocas.h"
#include "sparse_mat.h"
#include "ocas_helper.h"

#include "version.h"

void print_usage(void)
{
  printf("MSVMOCAS: Optimized Cutting Plane Algorithm for Support Vector Machines\n" 
         "          Training linear multi-class SVM classifier from examples\n"
         "          " OCAS_VERSION "\n"
         "\n"
         "   usage: msvmocas [options] example_file model_file\n"
         "\n"
         "Arguments:\n"
         "         example_file  file with training examples stored in SVM^light format\n"
         "         model_file    text file to store learned linear rule f(x)=W'*x; the model file\n"
         "                       contains M columns and D lines (M is number of classes, D is data\n"
         "                       dimension) corresponding to the elements of the matrix W [D x M].\n"
         "General options:\n"
         "         -h            this help\n"
         "         -v [0,1]      verbosity level (default 1).\n"
         "Learning options:\n"
         "         -c float      regularization constant C (default 1).\n"
         "         -n int        use only first n examples for training. By default n equals to\n"
         "                       the number of examples in the example_file.\n"
         "Optimization options:\n" 
         "         -m [0,1]      solver to be used: 0 ... standard cutting plane (BMRM, SVM^perf)\n"
         "                                          1 ... OCAS (default 1).\n"
         "         -s int        cache size for cutting planes (default 2000).\n"
         "Stopping conditions:\n"
         "         -a float      absolute tolerance TolAbs: halt if QP-QD <= TolAbs (default 0).\n"
         "         -r float      relative tolerance TolRel: halt if QP-QD <= abs(QP)*TolRel (default 0.01).\n"
         "         -q float      desired objective value QPValue: halt if QP <= QPValue (default 0).\n"         
         "         -t float      halts if the solver time (loading time is not counted) exceeds\n"
         "                       the given time in seconds (default inf).\n\n"
         "Example:\n"
         "  Train multi-class SVM classifier from example file ./data/example4_train.light with\n"
         "  regularization constant set to C = 10, verbosity switched off and save model to\n"
         "  file ./data/msvmocas.model\n"
         "    ./msvmocas -c 10 -v 0 ./data/example4_train.light ./data/msvmocas.model\n"
         "\n"  
         "  Compute testing error of the classifier stored in ./data/msvmocas.model using testing\n"
         "  examples from ./data/example4_test.light and save predicted labels to ./data/example4_test.pred\n"
         "    ./linclass -e -o ./data/exaple4_test.pred ./data/example4_test.light ./data/msvmocas.model\n"
         "\n"
         );
}


int main(int argc, char *argv[])
{
  double C, TolRel, TolAbs, QPBound, MaxTime;
  uint32_t i, j, BufSize;
  uint16_t Method;
  ocas_return_value_T ocas;
  int len;
  int recognized;
  int exitflag = 1;
  int verb;

  /* timing variables */
  double load_time;
  double total_time;

  char *model_fname;
  char *input_fname;
  FILE *fid;

  /* start time measuring */
  total_time = get_time();

  /* init */
  data_X = NULL;
  data_y = NULL;
  W = NULL;
  oldW = NULL;
  sparse_A.nz_dims = NULL;
  sparse_A.index = NULL;
  sparse_A.value = NULL;
  new_a = NULL;
  full_A = NULL;

  /* default setting of input arguments*/
  C = 1.0;
  Method = 1;
  TolRel = 0.01;
  TolAbs = 0.0;
  QPBound = 0.0;
  BufSize = 2000;
  MaxTime = (double)LIBOCAS_PLUS_INF;
  nData = -1;
  verb = 1;

  /*-----------------------------------------------------------
    Process input arguments 
  ------------------------------------------------------------*/
  if(argc ==1 || strcmp(argv[1], "-h") == 0)  
  {
    print_usage();
    goto clean_up;
  }

  if(argc < 3)
  {
    fprintf(stderr,"Not enough input arguments.\n\n");
    goto clean_up;
  }


  for (i = 1; i < argc-2; i++)  
  {
    recognized = 0;
    if (strcmp(argv[i], "-h") == 0)  
    {
      print_usage();
      goto clean_up;
    }

    if (strcmp(argv[i], "-c") == 0)  
    {
      if(i+1 >= argc-2)
      {
        fprintf(stderr,"You have to specify a value after argument -c\n"); 
        goto clean_up;  
      }
      C = atof(argv[i+1]); 
      if(C <=0)
      { 
        fprintf(stderr,"Parameter C must be geater than zero.\n"); 
        goto clean_up; 
      } 
      i++;
      recognized = 1;
    }

    if (strcmp(argv[i], "-n") == 0)  
    {
      if(i+1 >= argc-2)
      {
        fprintf(stderr,"You have to specify a value after argument -n\n"); 
        goto clean_up;  
      }
      nData = atol(argv[i+1]); 
      if(nData <=0)
      { 
        fprintf(stderr,"A value after the argument -n must be greater than zero.\n"); 
        goto clean_up; 
      } 
      i++;
      recognized = 1;
    }

    if (strcmp(argv[i], "-s") == 0)  
    {
      if(i+1 >= argc-2)
      {
        fprintf(stderr,"You have to specify a value after argument -s\n"); 
        goto clean_up;  
      }
      BufSize = atol(argv[i+1]); 
      if(nData <=0)
      { 
        fprintf(stderr,"A value after the argument -s must be greater than zero.\n"); 
        goto clean_up; 
      } 
      i++;
      recognized = 1;
    }


    if (strcmp(argv[i], "-m") == 0)  
    {
      if(i+1 >= argc-2)
      {
        fprintf(stderr,"You have to specify a value after argument -m\n"); 
        goto clean_up;  
      }
      Method = atoi(argv[i+1]); 
      if(Method != 0 && Method != 1)
      { 
        fprintf(stderr,"A value after the argument -m must be 0 or 1.\n"); 
        goto clean_up; 
      } 
      i++;
      recognized = 1;
    }

    if (strcmp(argv[i], "-v") == 0)  
    {
      if(i+1 >= argc-2)
      {
        fprintf(stderr,"You have to specify a value after argument -v\n"); 
        goto clean_up;  
      }
      verb = atoi(argv[i+1]); 
      if(verb < 0 || verb > 1)
      {
        fprintf(stderr,"A value after the argument -v must be either 0 or 1.\n"); 
        goto clean_up;  
      }
        

      i++;
      recognized = 1;
    }

    if (strcmp(argv[i], "-a") == 0)  
    {
      if(i+1 >= argc-2)
      {
        fprintf(stderr,"You have to specify a value after argument -a\n"); 
        goto clean_up;  
      }
      TolAbs = atof(argv[i+1]); 
      if(TolAbs < 0)
      { 
        fprintf(stderr,"A value after the argument -a must be a positive scalar.\n"); 
        goto clean_up; 
      } 
      i++;
      recognized = 1;
    }

    if (strcmp(argv[i], "-r") == 0)  
    {
      if(i+1 >= argc-2)
      {
        fprintf(stderr,"You have to specify a value after argument -r\n"); 
        goto clean_up;  
      }
      TolRel = atof(argv[i+1]); 
      if(TolRel < 0)
      { 
        fprintf(stderr,"A value after the argument -r must be a positive scalar.\n"); 
        goto clean_up; 
      } 
      i++;
      recognized = 1;
    }

    if (strcmp(argv[i], "-q") == 0)  
    {
      if(i+1 >= argc-2)
      {
        fprintf(stderr,"You have to specify a value after argument -q\n"); 
        goto clean_up;  
      }
      QPBound = atof(argv[i+1]); 
      i++;
      recognized = 1;
    }

    if (strcmp(argv[i], "-t") == 0)  
    {
      if(i+1 >= argc-2)
      {
        fprintf(stderr,"You have to specify a value after argument -t\n"); 
        goto clean_up;  
      }
      MaxTime = atof(argv[i+1]); 
      if(MaxTime <=0)
      {
        fprintf(stderr,"A value after the argument -t must be a positive scalar.\n"); 
        goto clean_up; 
      } 
      i++;
      recognized = 1;
    }

    if(recognized == 0)
    {
      fprintf(stderr,"Unknown input argument: %s\n", argv[i]);
      goto clean_up;  
    }
        
  }

  len = strlen(argv[argc-2]);
  input_fname = calloc(len+1,sizeof(char));
  strcpy(input_fname, argv[argc-2]);

  len = strlen(argv[argc-1]);
  model_fname = calloc(len+1,sizeof(char));
  strcpy(model_fname, argv[argc-1]);

  if(verb)
    printf("Input file: %s\n", input_fname);

  /*----------------------------------------------------------------
    Load input examples
  -------------------------------------------------------------------*/
  load_time = get_time();
  if( load_svmlight_file(input_fname,verb) == -1 || data_X == NULL || data_y == NULL)
    goto clean_up;

  load_time = get_time() - load_time;
  
  /* get examples' dimension */
  nDim = mxGetM(data_X);

  /* if not given set number of training examples to be used */
  if(nData == -1)
    nData = mxGetN(data_X);
  else if(nData <= 0 || nData > mxGetN(data_X))
  {
    fprintf(stderr,"Number of examples in the input file is %d.\n", mxGetN(data_X));
    fprintf(stderr,"A value after argument -n must be less or equal to the number of examples.\n");
    goto clean_up;
  }

  /* get number of classes */
  for(i=0, nY = 0; i < nData; i++)
    nY = LIBOCAS_MAX(nY, (uint32_t)data_y[i]);


  /*----------------------------------------------------------------
    Print setting
  -------------------------------------------------------------------*/
  if(verb)
  {
    printf("Input file statistics:\n"
           "   # of examples  : %d\n"
           "   # of classes   : %d\n"
           "   dimensionality : %d\n",
           mxGetN(data_X), nY, nDim);
    if( mxIsSparse(data_X) )
      printf("   density        : %.2f%% (sparse matrix representation used)\n",
             100.0*(double)mxGetNZMAX(data_X)/((double)nDim*(double)(mxGetN(data_X))));
    else
      printf("   density        : 100%% (full matrix representation used)\n");
    
    printf("Setting:\n"
           "   C              : %f\n"
           "   # of examples  : %d\n"
           "   solver         : %d\n"
           "   cache size     : %d\n"
           "   TolAbs         : %f\n"
           "   TolRel         : %f\n"
           "   QPValue        : %f\n"
           "   MaxTime        : %f [s]\n"
           "   Verbosity      : %d\n",
           C, nData, Method,BufSize,TolAbs,TolRel, QPBound, MaxTime,verb);
  }


  /*----------------------------------------------------------------
    Allocate memory for working variables and cutting plane cache
  -------------------------------------------------------------------*/

  /* learned weight vector */
  W = (double*)mxCalloc(nDim*nY,sizeof(double));
  if(W == NULL)
  {
    fprintf(stderr,"Not enough memory for matrix W.\n");
    goto clean_up;
  }
    
  oldW = (double*)mxCalloc(nDim*nY,sizeof(double));
  if(oldW == NULL) 
  {
    fprintf(stderr,"Not enough memory for matrix oldW.");
    goto clean_up;
  }

  /* allocate buffer for computing cutting plane */
  new_a = (double*)mxCalloc(nDim*nY,sizeof(double));
  if(new_a == NULL) 
  {
    fprintf(stderr,"Not enough memory for auxciliary cutting plane buffer new_a.");  
    goto clean_up;
  }

  if(mxIsSparse(data_X))
  {

    /* init cutting plane buffer */
    sparse_A.nz_dims = mxCalloc(BufSize,sizeof(uint32_t));
    sparse_A.index = mxCalloc(BufSize,sizeof(sparse_A.index[0]));
    sparse_A.value = mxCalloc(BufSize,sizeof(sparse_A.value[0]));
    if(sparse_A.nz_dims == NULL || sparse_A.index == NULL || sparse_A.value == NULL) 
    {
        fprintf(stderr,"Not enough memory for cutting plane buffer sparse_A.");  
        goto clean_up;
    }

    /*----------------------------------------------------------------
      Run OCAS run ...
      -------------------------------------------------------------------*/
    if(verb)
    {
      printf("Starting optimization:\n");

      ocas = msvm_ocas_solver( C, data_y, nY, nData, TolRel, TolAbs, QPBound, MaxTime,BufSize, Method,
                               &msvm_sparse_compute_W, &msvm_full_update_W, &msvm_sparse_add_new_cut,
                               &msvm_sparse_compute_output, &qsort_data, &ocas_print, 0);
    }
    else
    {
      ocas = msvm_ocas_solver( C, data_y, nY, nData, TolRel, TolAbs, QPBound, MaxTime,BufSize, Method,
                               &msvm_sparse_compute_W, &msvm_full_update_W, &msvm_sparse_add_new_cut,
                               &msvm_sparse_compute_output, &qsort_data, &ocas_print_null, 0);
    }
  }
  else
  {
    /* init cutting plane buffer */
    full_A = mxCalloc(BufSize*nDim*nY,sizeof(double));
    if( full_A == NULL )
    {
      fprintf(stderr,"Not enough memory for cutting plane buffer full_A.");  
      goto clean_up;
    }

    if(verb)
    {
      printf("Starting optimization:\n");
    
      ocas = msvm_ocas_solver( C, data_y, nY, nData, TolRel, TolAbs, QPBound, MaxTime,BufSize, Method,
                               &msvm_full_compute_W, &msvm_full_update_W, &msvm_full_add_new_cut,
                               &msvm_full_compute_output, &qsort_data, &ocas_print, 0); 
    }
    else
    {
      ocas = msvm_ocas_solver( C, data_y, nY, nData, TolRel, TolAbs, QPBound, MaxTime,BufSize, Method,
                               &msvm_full_compute_W, &msvm_full_update_W, &msvm_full_add_new_cut,
                               &msvm_full_compute_output, &qsort_data, &ocas_print_null, 0); 

    }

  }


  total_time=get_time()-total_time;

  if(verb)
  {
    printf("Stopping condition: ");
    switch( ocas.exitflag )
    {
       case 1: printf("1-Q_D/Q_P <= TolRel(=%f) satisfied.\n", TolRel); break;
       case 2: printf("Q_P-Q_D <= TolAbs(=%f) satisfied.\n", TolAbs); break;
       case 3: printf("Q_P <= QPBound(=%f) satisfied.\n", QPBound); break;
       case 4: printf("Optimization time (=%f) >= MaxTime(=%f).\n", ocas.ocas_time, MaxTime); break;
       case -1: printf("Has not converged!\n" ); break;
       case -2: printf("Not enough memory for the solver.\n" ); break;
    }

    printf("Timing statistics:\n"
           "   load_time      : %f[s]\n"
           "   qp_solver_time : %f[s]\n"
           "   sort_time      : %f[s]\n"
           "   output_time    : %f[s]\n"
           "   add_time       : %f[s]\n"
           "   w_time         : %f[s]\n"
           "   print_time     : %f[s]\n"
           "   ocas_time      : %f[s]\n"
           "   total_time     : %f[s]\n",
           load_time, ocas.qp_solver_time, ocas.sort_time, ocas.output_time, 
           ocas.add_time, ocas.w_time, ocas.print_time, ocas.ocas_time, total_time);

    printf("Training error: %.4f%%\n", 100*(double)ocas.trn_err/(double)nData);
  }

  /* save W to the model file */
  fid = fopen(model_fname, "w+");
  if(fid == NULL) {
    perror("fopen error ");
    mexPrintf("Cannot open model file.");
    goto clean_up;
  }

  for(i=0; i < nDim; i++) 
  {
    for(j=0; j < nY; j++)
    {
      fprintf(fid, "%.20f ", W[LIBOCAS_INDEX(i,j,nDim)]);
    }
    fprintf(fid, "\n");
  }

  fclose(fid);

  if(verb)
    printf("Result saved to model file: %s\n", model_fname);

  exitflag = 0;

clean_up:

  mxDestroyArray(data_X);
  mxFree(data_y);
  mxFree(W);
  mxFree(oldW);
  mxFree(new_a);
  mxFree(full_A);

  mxFree(sparse_A.nz_dims);
  if( sparse_A.index !=NULL) 
  {
    for(i=0; i < BufSize; i++)
      if(sparse_A.index[i] != NULL)
        mxFree(sparse_A.index[i]);

    mxFree(sparse_A.index);
  }
  if( sparse_A.value != NULL)
  {
    for(i=0; i < BufSize; i++)
      if(sparse_A.value[i] != NULL)
        mxFree(sparse_A.value[i]);
    
    mxFree(sparse_A.value);
  }

  return(exitflag);
}



