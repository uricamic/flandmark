/*-----------------------------------------------------------------------
 * linclass.c: Implementation of linear classification rule classifying 
 *  examples from the SVM^light format.
 *   
 * Copyright (C) 2008, 2009 Vojtech Franc, xfrancv@cmp.felk.cvut.cz
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
#include <stdint.h>

#include "lib_svmlight_format.h"
#include "libocas.h"
#include "version.h"

#define MODELFILE_MAXLINELEN 1000000

void print_usage(void)
{
  printf("LINCLASS: Predict labels by linear classication rule\n" 
         "          " OCAS_VERSION "\n"
         "\n"
         "   usage: linclass [options] example_file model_file\n"
         "\n"
         "Arguments:\n"
         "         example_file    text file with testing examples stored in SVM^light format\n"
         "         model_file      text file which contains either binary (two-class) linear\n"
         "                         rule f(x)=w'*x+w0 or multi-class rule f(x)=W'*x which are\n"
         "                         produced by svmocas and msvmocas, respectively\n"
         "Options:\n"
         "         -e              prints classification error computed from predicted\n"
         "                         labels and labels contained in the example_file.\n"
         "         -h              this help\n"
         "         -o output_file  save predictions to output_file rather than to stdout.\n"
         "         -t [0,1]        output type: 0 .. predicted labels (default)\n"
         "                                      1 .. discriminant values\n"
         "         -v [0,1]        verbosity level (default 0).\n"
         "\n"
         "Examples\n"
         "  Train SVM classifier from riply_trn.light with regularization constant C = 10,\n"
         "  bias switched on, verbosity switched off and model saved to svmocas.model\n"
         "    ./svmocas -c 10 -b 1 -v 0 riply_trn.light svmocas.model \n"
         "\n"
         "  Compute testing error of the classifier stored in svmocas.model using testing\n"
         "  examples from riply_tst.light and save predicted labels to riply_tst.pred\n"
         "    ./linclass -e -o riply_tst.pred riply_tst.light svmocas.model\n"
         "\n"
         );
}

int main(int argc, char *argv[])
{
  uint32_t i, j;
  int len;
  int recognized;
  int exitflag = 0;
  int verb;
  int binary_problem;
  int output_type;
  int print_error;
  char *line;

  char *model_fname;
  char *input_fname;
  char *output_fname;
  FILE *fid, *fout;
  
  double *feat_val;
  uint32_t *feat_idx;

  int go = 1;
  char *endptr, *begptr;
  int nLines = 0;
  int nCols = 0, tmp_nCols;
  double val;
  double *W;
  double W0;
  uint32_t nY, nDim;

  long line_cnt = 0;
  int label, pred_label=0;
  long max_dim = 0;
  long nnzf;
  double dfce, max_dfce;
  long nErrors = 0;
  long *nClassErrors, *nClass;

  /* init */
  fid = NULL;
  fout = NULL;
  feat_val = NULL;
  feat_idx = NULL;
  W = NULL;
  W0 = 0;
  output_fname = NULL;
  nClassErrors = NULL;
  nClass = NULL;

  /* default setting of input arguments*/
  verb = 0;
  print_error = 0;
  output_type = 0;  

  /* Allocate memory */
  line = calloc(MODELFILE_MAXLINELEN, sizeof(char));
  if( line == NULL )
  {
    fprintf(stderr,"Not enough memmory to allocate line buffer.\n");
    goto clean_up;
  }

  feat_idx = calloc(LIBSLF_MAXLINELEN, sizeof(uint32_t));
  if( feat_idx == NULL )
  {
    fprintf(stderr,"Not enough memmory to allocate feat_idx.\n");
    goto clean_up;
  }

  feat_val = calloc(LIBSLF_MAXLINELEN, sizeof(double));
  if( feat_val == NULL )
  {
    fprintf(stderr,"Not enough memmory to allocate feat_val.\n");
    goto clean_up;
  }
  

  /*-----------------------------------------------------------
    Process input arguments 
  ------------------------------------------------------------*/
  if(argc ==1 || strcmp(argv[1], "-h") == 0)  
  {
    print_usage();
    goto clean_up;
  }

  if(argc < 2)
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

    if (strcmp(argv[i], "-e") == 0)  
    {
      print_error = 1;
      recognized = 1;
      continue;
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
      continue;
    }

    if (strcmp(argv[i], "-t") == 0)  
    {
      if(i+1 >= argc-2)
      {
        fprintf(stderr,"You have to specify a value after argument -t\n"); 
        goto clean_up;  
      }
      output_type = atoi(argv[i+1]); 
      if(output_type != 0 && output_type != 1)
      {
        fprintf(stderr,"A value after the argument -t must be either 0 or 1.\n"); 
        goto clean_up;  
      }
        
      i++;
      recognized = 1;
      continue;
    }

    if (strcmp(argv[i], "-o") == 0)  
    {
      if(i+1 >= argc-2)
      {
        fprintf(stderr,"You have to specify a string after argument -o\n"); 
        goto clean_up;  
      }

      len = strlen(argv[i+1]);
      output_fname = calloc(len+1,sizeof(char));
      strcpy(output_fname, argv[i+1]);

      i++;
      recognized = 1;
      continue;
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
  {
    printf("Verbosity: %d\n", verb);
    printf("Output type: %d\n", output_type);
    printf("Print error: %d\n", print_error);
    printf("Example file: %s\n", input_fname);
    printf("Model file: %s\n", model_fname);
    if( output_fname != NULL)
      printf("Output file: %s\n", output_fname);
    else
      printf("Output file: stdout\n");
  }


  /*----------------------------------------------------------------
    Load classification rule which is either 
     vector [nDim x 1] + bias [1x1]
    or 
     matrix [nDim x nY]
  -------------------------------------------------------------------*/
  
  /* load W from model file */
  fid = fopen(model_fname, "r");
  if(fid == NULL) {
    fprintf(stderr,"Cannot open model file.\n");
    perror("fopen error ");
    goto clean_up;
  }
  
  if(verb)
  {
    printf("Analysing model file... ");
    fflush(stdout);
  }

  /* read the first line */
  if(fgets(line,LIBSLF_MAXLINELEN, fid) == NULL ) 
  {
    fprintf(stderr,"Empty example file.\n");
    goto clean_up;
  }
  else
  {
    nLines = 1;
    begptr = line;
    while(1)
    {
      val = strtod(begptr, &endptr);

      if(val == 0 && begptr == endptr)
        break;

      nCols++;
      begptr = endptr;
    }
  }

  go = 1;
  while(go) 
  {
    begptr = line;

    tmp_nCols = 0;
    while(1)
    {
      val = strtod(begptr, &endptr);

      if(val == 0 && begptr == endptr)
        break;

      tmp_nCols++;
      begptr = endptr;
    }
    if( tmp_nCols != nCols)
    {
      fprintf(stderr,"Error: Model file contains lines with different number of colums.\n");
      goto clean_up;
    }

    if(fgets(line,LIBSLF_MAXLINELEN, fid) == NULL ) 
    {
      go = 0;
    }
    else
      nLines++;
  }
  
  if(verb)
    printf("done.\n"
           "Number of lines: %d\n"
           "Number of columns: %d\n",
           nLines,nCols);

  if(nCols == 1)
  {
    nY = 2;
    nDim = nLines-1;
    binary_problem = 1;

    /* learned weight vector */
    W = (double*)calloc(nDim,sizeof(double));
    if(W == NULL)
    {
      fprintf(stderr,"Not enough memory for vector W.\n");
      goto clean_up;
    }

    if(verb)
    {
      printf("Model file contains binary classification rule.\n");
      printf("Reading model file...");
    }

    fseek(fid,0,SEEK_SET);
    for(i=0; i <= nDim; i++)
    {
      if(fgets(line,LIBSLF_MAXLINELEN, fid) == NULL ) 
      {
        fprintf(stderr,"Model file corrupted.\n");
        goto clean_up;
      }

      begptr = line;
      val = strtod(begptr, &endptr);

      if(val == 0 && begptr == endptr)
      {
        fprintf(stderr,"Model file corrupted.\n");
        goto clean_up;
      }

      if(i < nDim)
        W[i] = val;
      else
        W0 = val;
    }

    if(verb)
      printf("done.\n");

  }
  else
  {
    nY = nCols;
    nDim = nLines;
    binary_problem = 0;

    /* learned weight vector */
    W = (double*)calloc(nDim*nY,sizeof(double));
    if(W == NULL)
    {
      fprintf(stderr,"Not enough memory for matrix W.\n");
      goto clean_up;
    }

    if(verb)
    {
      printf("Model file contains multi-class classification rule.\n");
      printf("Reading model file...");
    }

    fseek(fid,0,SEEK_SET);
    for(i=0; i < nDim; i++)
    {
      if(fgets(line,LIBSLF_MAXLINELEN, fid) == NULL ) 
      {
        fprintf(stderr,"Model file corrupted.\n");
        goto clean_up;
      }

      begptr = line;
      for(j=0; j < nY; j++)
      {
        val = strtod(begptr, &endptr);

        if(val == 0 && begptr == endptr)
        {
          fprintf(stderr,"Model file corrupted.\n");
          goto clean_up;
        }
        begptr = endptr;
        
        
        W[LIBOCAS_INDEX(i,j,nDim)] = val;
      }
    }
    if(verb)
      printf("done.\n");
  }
  
  fclose(fid);

/*  printf("W0=%f, W = [ ", W0);*/
/*  for(i=0; i < nDim; i++)*/
/*    printf("%f ", W[i]);*/
/*  printf("]\n");*/
/*  printf("W = [\n");*/
/*  for(j=0; j < nDim; j++)*/
/*  {*/
/*    for(i=0; i < nY; i++)*/
/*      printf("%f ", W[LIBOCAS_INDEX(j,i,nDim)]);*/
/*    printf("\n");*/
/*  }*/

/*  load_time = get_time() - load_time;*/

  /*-----------------------------------------------------
    Read examples and classify them.
    -----------------------------------------------------*/

  fid = fopen(input_fname, "r");
  if(fid == NULL) {
    fprintf(stderr,"Cannot open input file.\n");
    perror("fopen error ");
    goto clean_up;
  }

  if(output_fname == NULL)
    fout = stdout;
  else
  {
    fout = fopen(output_fname, "w+");
    if(fid == NULL) {
      fprintf(stderr,"Cannot open output file.\n");
      perror("fopen error ");
      fclose(fid);
      goto clean_up;
    }
  }

  if(verb)
  { 
    if(output_fname != NULL)
      printf("Classifying...");
    else
      printf("Outputs:\n");
  }

  nClassErrors = (long*)calloc(nY,sizeof(long));
  if(nClassErrors == NULL)
  {
    fprintf(stderr,"Not enough memory for vector nClassError.\n");
    goto clean_up;
  }

  nClass = (long*)calloc(nY,sizeof(long));
  if(nClass == NULL)
  {
    fprintf(stderr,"Not enough memory for vector nClass.\n");
    goto clean_up;
  }
  

  go = 1;
  while(go) { 
    
    if(fgets(line,LIBSLF_MAXLINELEN, fid) == NULL ) 
    {
      go = 0;
    }
    else
    {
      line_cnt ++;
      nnzf = svmlight_format_parse_line(line, &label, feat_idx, feat_val);
      
      if(nnzf == -1) 
      {
         fprintf(stderr,"Parsing error on line %ld .\n", line_cnt);
         fprintf(stderr,"Probably defective input file.\n");
         goto clean_up;
      }

      max_dim = LIBOCAS_MAX(max_dim,feat_idx[nnzf-1]);

      if(binary_problem == 1)
      {
        dfce = W0;
        for(i=0; i < nnzf; i++)
        {
          if(feat_idx[i]-1 < nDim)
            dfce += feat_val[i]*W[feat_idx[i]-1];
        }

        if(label == +1)
          nClass[0]++;
        else
          nClass[1]++;

        if(dfce >=0 && label == -1)
        {
          nClassErrors[1]++;
          nErrors++;
        }
        else if (dfce < 0 && label== +1)
        {
          nClassErrors[0]++;
          nErrors++;
        }

        if(output_type == 0)
        {
          if(dfce >=0 )
            fprintf(fout,"+1\n");
          else
            fprintf(fout,"-1\n");
            
        }
        else
          fprintf(fout,"%.20f\n", dfce);
      }
      else
      {
        max_dfce = -LIBOCAS_PLUS_INF;
        for(j=0; j < nY; j++)
        {
          dfce = 0;
          for(i=0; i < nnzf; i++)
          {
            if(feat_idx[i]-1 < nDim)
              dfce += feat_val[i]*W[LIBOCAS_INDEX(feat_idx[i]-1,j,nDim)];
          }
          if(output_type==1)
            fprintf(fout,"%.20f ", dfce);

          if(max_dfce < dfce)
          {
            max_dfce = dfce;
            pred_label = j+1;
          }
        }
        if(output_type==0)
            fprintf(fout,"%d", pred_label);

        fprintf(fout,"\n");

        nClass[label-1]++;

        if(label != pred_label)
        {
          nErrors++;
          nClassErrors[label-1]++;
        }

      }
    }
  }  

  if(verb)
  {
    if(output_fname != NULL)
      printf("done.\n");

    printf("Number of examples: %ld\n"
           "Maximal dimensionality: %ld\n", line_cnt, max_dim);
  }
  if(print_error)
  {
    printf("Classification error: %f%%(%ld/%ld)\n", 100.0*(double)nErrors/(double)line_cnt,nErrors,line_cnt);
    printf("Per-class errors: ");
    if(binary_problem)
    {
      printf("+1: %f%%(%ld/%ld) -1: %f%%(%ld/%ld)\n", 
             100.0*(double)nClassErrors[0]/(double)nClass[0], nClassErrors[0],nClass[0],
             100.0*(double)nClassErrors[1]/(double)nClass[1], nClassErrors[1],nClass[1]);
    }
    else 
    {
      for(i=0; i < nY; i++)
        printf("%d: %f%%(%ld/%ld) ", i+1, 100.0*(double)nClassErrors[i]/(double)nClass[i], 
               nClassErrors[i],nClass[i]);
      printf("\n");
    } 
  } 


  fclose(fid);
  fclose(fout);

  exitflag = 1;
  
clean_up:

  free(W);
  free(line);
  free(feat_val);
  free(feat_idx);
  free(nClassErrors);
  free(nClass);

  return(exitflag);
}



