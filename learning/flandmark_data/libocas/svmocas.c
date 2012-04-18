/*-----------------------------------------------------------------------
 * svmocas.c: Standalone application implementing the OCAS solver for 
 *   training linear SVM classifiers.
 *   
 * Copyright (C) 2008, 2009 Vojtech Franc, xfrancv@cmp.felk.cvut.cz
 *                    Soeren Sonnenburg, soeren.sonnenburg@first.fraunhofer.de
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public 
 * License as published by the Free Software Foundation; 
 * Version 3, 29 June 2007
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
  printf("SVMOCAS: Optimized Cutting Plane Algorithm for Support Vector Machines\n"
         "         Training binary linear SVM classifier from examples\n"
         "         " OCAS_VERSION "\n"
         "\n"
         "   usage: svmocas [options] example_file model_file\n"
         "\n"
         "Arguments:\n"
         "         example_file  file with training examples stored in SVM^light format\n"
         "         model_file    file to store learned linear rule f(x)=w'*x+w0; the model file\n"
         "                       contains d+1 lines (d is data dimension); first n lines are\n"
         "                       coordinates of w and the last one is w0.\n"
         "General options:\n"
         "         -h            this help\n"
         "         -v [0,1]      verbosity level (default 1).\n"
         "Learning options:\n"
         "         -c float      regularization constant C (default 1).\n"
         "         -C filename   if specified each example has different reg. constant given in the text file.\n"
         "                       Each line of the text file must contain a single constant (positive double)\n"
         "                       for the corresponding example. If -C is used then the option -c is ignored.\n"
         "         -b [0,1]      value of L2-bias feature. A value of 0 implies not having bias (default 0).\n"
         "         -n int        use only first n examples for training. By default n equals to\n"
         "                       the number of examples in the example_file.\n"
         "Optimization options:\n" 
         "         -m [0,1]      solver to be used: 0 ... standard cutting plane (equivalent to BMRM, SVM^perf)\n"
         "                                          1 ... OCAS (default 1).\n"
         "         -s int        cache size for cutting planes (default 2000).\n"
         "         -p int        number of threads (default 1).\n"
         "Stopping conditions:\n"
         "         -a float      absolute tolerance TolAbs: halt if QP-QD <= TolAbs (default 0).\n"
         "         -r float      relative tolerance TolRel: halt if QP-QD <= abs(QP)*TolRel (default 0.01).\n"
         "         -q float      desired objective value QPValue: halt is QP <= QPValue (default 0).\n"         
         "         -t float      halts if the solver time (loading time is not counted) exceeds\n"
         "                       the given time in seconds (default inf).\n\n"
         "Example:\n"
         "  Train binary SVM classifier from ./data/riply_trn.light with regularization constant C = 10,\n"
         "  bias switched on, verbosity switched off and save model to ./data/svmocas.model\n"
         "    ./svmocas -c 10 -b 1 -v 0 ./data/riply_trn.light ./data/svmocas.model \n"
         "\n"
         "  Compute testing error of the classifier stored in ./data/svmocas.model using testing\n"
         "  examples from ./data/riply_tst.light and save predicted labels to ./data/riply_tst.pred\n"
         "    ./linclass -e -o ./data/riply_tst.pred ./data/riply_tst.light ./data/svmocas.model\n"
         "\n"
         );
}


int main(int argc, char *argv[])
{
  double C, TolRel, TolAbs, QPBound, MaxTime;
  uint32_t i, j, BufSize;
  uint16_t Method;
  int number_of_threads;
  ocas_return_value_T ocas;
  int len;
  int recognized;
  int exitflag = 1;
  int verb;

  double *vec_C = NULL;
  uint32_t len_vec_C =0;

  /* timing variables */
  double load_time;
  double total_time;

  char *model_fname = NULL;
  char *input_fname = NULL;
  char *regconst_fname = NULL;
  FILE *fid;

  /* start time measuring */
  total_time = get_time();

  /* init */
  data_X = NULL;
  data_y = NULL;
  W = NULL;
  oldW = NULL;
  A0 = NULL;
  sparse_A.nz_dims = NULL;
  sparse_A.index = NULL;
  sparse_A.value = NULL;
  new_a = NULL;
  full_A = NULL;

  /* default setting of input arguments*/
  X0 = 0;
  C = 1.0;
  Method = 1;
  TolRel = 0.01;
  TolAbs = 0.0;
  QPBound = 0.0;
  BufSize = 2000;
  MaxTime = (double)LIBOCAS_PLUS_INF;
  nData = -1;
  number_of_threads = 1;
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

    if (strcmp(argv[i], "-C") == 0)
    {
      if(i+1 >= argc-2)
      {
        fprintf(stderr,"You have to specify a file name after argument -C\n"); 
        goto clean_up;  
      }

      len = strlen( argv[i+1] );
      regconst_fname = calloc(len+1,sizeof(char));
      strcpy(regconst_fname, argv[i+1]);
      
      i++;
      recognized = 1;
    }

    if (strcmp(argv[i], "-p") == 0)  
    {
      if(i+1 >= argc-2)
      {
        fprintf(stderr,"You have to specify a value after argument -p\n"); 
        goto clean_up;  
      }
      number_of_threads = atoi(argv[i+1]); 
      if(number_of_threads <= 0)
      { 
        fprintf(stderr,"Argumnet after parameter -p must be geater than zero.\n"); 
        goto clean_up; 
      } 
      i++;
      recognized = 1;
    }


    if (strcmp(argv[i], "-b") == 0)  
    {
      if(i+1 >= argc-2)
      {
        fprintf(stderr,"You have to specify a value after argument -b\n"); 
        goto clean_up;  
      }
      X0 = atof(argv[i+1]); 
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

  /* load regularization constants from file if requested */
  if( regconst_fname != NULL)
  {
    if(load_regconsts(regconst_fname, &vec_C, &len_vec_C, verb) == -1)
    {
      fprintf(stderr,"Cannot load regularization constants from file %s\n", regconst_fname); 
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


  /*----------------------------------------------------------------
    Print setting
  -------------------------------------------------------------------*/
  if(verb)
  {
    printf("Input file statistics:\n"
           "   # of examples  : %d\n"
           "   dimensionality : %d\n",
           mxGetN(data_X), nDim);
    if( mxIsSparse(data_X) )
      printf("   density        : %.2f%% (sparse matrix representation used)\n",
             100.0*(double)mxGetNZMAX(data_X)/((double)nDim*(double)(mxGetN(data_X))));
    else
      printf("   density        : 100%% (full matrix representation used)\n");
     
    printf("Setting:\n");
             
    if(vec_C == NULL)
      printf("   C              : %f\n", C);
    else
    {
      printf("   C              : different for each example; specified in %s\n", regconst_fname);
      for(i = 0; i < len_vec_C; i++)
        printf(" %f", vec_C[i]);
      printf("\n");
    }

    printf("   bias           : %.0f\n"
           "   # of examples  : %d\n"
           "   solver         : %d\n"
           "   cache size     : %d\n"
           "   # of threads   : %d\n"
           "   TolAbs         : %f\n"
           "   TolRel         : %f\n"
           "   QPValue        : %f\n"
           "   MaxTime        : %f [s]\n"
           "   verbosity      : %d\n",
           X0, nData, Method,BufSize,number_of_threads,TolAbs,TolRel, QPBound, MaxTime,verb);
  }


  if( vec_C != NULL && len_vec_C < nData)
  {
    fprintf(stderr, "Number of regularization constants (%d) is less then number of examples (%d).\n",
            len_vec_C, nData);
    goto clean_up;
  }

  /*----------------------------------------------------------------
    Allocate memory for working variables and cutting plane cache
  -------------------------------------------------------------------*/

  /* learned weight vector */
  W = (double*)mxCalloc(nDim,sizeof(double));
  if(W == NULL)
  {
    fprintf(stderr,"Not enough memory for vector W.\n");
    goto clean_up;
  }
    
  oldW = (double*)mxCalloc(nDim,sizeof(double));
  if(oldW == NULL) 
  {
    fprintf(stderr,"Not enough memory for vector oldW.");
    goto clean_up;
  }

  W0 = 0;
  oldW0 = 0;

  A0 = mxCalloc(BufSize,sizeof(A0[0]));
  if(A0 == NULL) 
  {
    fprintf(stderr,"Not enough memory for vector A0.");
    goto clean_up;
  }

  /* allocate buffer for computing cutting plane */
  new_a = (double*)mxCalloc(nDim*number_of_threads,sizeof(double));
  if(new_a == NULL) 
  {
    fprintf(stderr,"Not enough memory for auxciliary cutting plane buffer new_a.");  
    goto clean_up;
  }


  if(vec_C != NULL)
  {
    for(i=0; i < nData; i++) 
      data_y[i] = data_y[i]*vec_C[i];
  }


  if(mxIsSparse(data_X))
  {

    /* for i=1:nData, X(:,i) = X(:,i)*y(i); end*/
    for(i=0; i < nData; i++) 
      mul_sparse_col(data_y[i], data_X, i);

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
    if(number_of_threads == 1)
    {
      if(verb)
      {
        printf("Starting optimization:\n");
        if( vec_C == NULL)
          ocas = svm_ocas_solver( C, nData, TolRel, TolAbs, QPBound, MaxTime,BufSize, Method, 
                                  &sparse_compute_W, &sparse_update_W, &sparse_add_new_cut, 
                                  &sparse_compute_output, &qsort_data, &ocas_print, 0);
        else
          ocas = svm_ocas_solver_difC( vec_C, nData, TolRel, TolAbs, QPBound, MaxTime,BufSize, Method, 
                                       &sparse_compute_W, &sparse_update_W, &sparse_add_new_cut, 
                                       &sparse_compute_output, &qsort_data, &ocas_print, 0);

      }
      else
      {
        if( vec_C == NULL)          
          ocas = svm_ocas_solver( C, nData, TolRel, TolAbs, QPBound, MaxTime,BufSize, Method, 
                                  &sparse_compute_W, &sparse_update_W, &sparse_add_new_cut, 
                                  &sparse_compute_output, &qsort_data, &ocas_print_null, 0);
        else
          ocas = svm_ocas_solver_difC( vec_C, nData, TolRel, TolAbs, QPBound, MaxTime,BufSize, Method, 
                                       &sparse_compute_W, &sparse_update_W, &sparse_add_new_cut, 
                                       &sparse_compute_output, &qsort_data, &ocas_print_null, 0);
          
      }
    }
    else
    {

      if(init_parallel_ocas(number_of_threads) != 0)
      {
        goto clean_up;
      }

      if(verb)
      {
        printf("Starting optimization:\n");
    
        if( vec_C == NULL)          
          ocas = svm_ocas_solver( C, nData, TolRel, TolAbs, QPBound, MaxTime,BufSize, Method, 
                                  &sparse_compute_W, &sparse_update_W, &parallel_sparse_add_new_cut, 
                                  &parallel_sparse_compute_output, &parallel_qsort_data, &ocas_print, 0);
        else
          ocas = svm_ocas_solver_difC( vec_C, nData, TolRel, TolAbs, QPBound, MaxTime,BufSize, Method, 
                                       &sparse_compute_W, &sparse_update_W, &parallel_sparse_add_new_cut, 
                                       &parallel_sparse_compute_output, &parallel_qsort_data, &ocas_print, 0);

      }
      else
      {
        if( vec_C == NULL)          
          ocas = svm_ocas_solver( C, nData, TolRel, TolAbs, QPBound, MaxTime,BufSize, Method, 
                                  &sparse_compute_W, &sparse_update_W, &parallel_sparse_add_new_cut, 
                                  &parallel_sparse_compute_output, &parallel_qsort_data, &ocas_print_null, 0);
        else
          ocas = svm_ocas_solver_difC( vec_C, nData, TolRel, TolAbs, QPBound, MaxTime,BufSize, Method, 
                                  &sparse_compute_W, &sparse_update_W, &parallel_sparse_add_new_cut, 
                                  &parallel_sparse_compute_output, &parallel_qsort_data, &ocas_print_null, 0);
          
      }
      destroy_parallel_ocas();
    }
  }
  else
  {
    if(number_of_threads != 1 && verb)
    {
      printf("Warning: Multi-thread support works currently only for sparse data.\n");
    }

    double *ptr = mxGetPr(data_X);
    for(i=0; i < nData; i++) {
      for(j=0; j < nDim; j++ ) {
        ptr[LIBOCAS_INDEX(j,i,nDim)] = ptr[LIBOCAS_INDEX(j,i,nDim)]*data_y[i];
      }
    }

    /* init cutting plane buffer */
    full_A = mxCalloc(BufSize*nDim,sizeof(double));
    if( full_A == NULL )
    {
      fprintf(stderr,"Not enough memory for cutting plane buffer full_A.");  
      goto clean_up;
    }

    if(verb)
    {
      printf("Starting optimization:\n");

      if(vec_C == NULL)
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
      if(vec_C == NULL)
        ocas = svm_ocas_solver( C, nData, TolRel, TolAbs, QPBound, MaxTime,BufSize, Method, 
                                &full_compute_W, &full_update_W, &full_add_new_cut, 
                                &full_compute_output, &qsort_data, &ocas_print_null, 0);
      else
        ocas = svm_ocas_solver_difC( vec_C, nData, TolRel, TolAbs, QPBound, MaxTime,BufSize, Method, 
                                     &full_compute_W, &full_update_W, &full_add_new_cut, 
                                     &full_compute_output, &qsort_data, &ocas_print_null, 0);

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

  /* save W and W0 to the model file */
  fid = fopen(model_fname, "w+");
  if(fid == NULL) {
    perror("fopen error ");
    mexPrintf("Cannot open model file.");
    goto clean_up;
  }

  for(i=0; i < nDim; i++)
    fprintf(fid, "%.20f\n", W[i]);

  fprintf(fid, "%.20f\n", W0);
  fclose(fid);

  if(verb)
    printf("Result saved to model file: %s\n", model_fname);

  exitflag = 0;

clean_up:

  mxDestroyArray(data_X);
  mxFree(data_y);
  mxFree(W);
  mxFree(oldW);
  mxFree(A0);
  mxFree(new_a);
  mxFree(full_A);
  mxFree(regconst_fname);
  mxFree(model_fname);
  mxFree(input_fname);
  mxFree(vec_C);

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



