/*=================================================================
 * SVMLIGHT_LINCLASS classifies examples in SVM^light file by linear rule.
 *
 *  Synopsis:
 *    [score,true_lab] = svmlight_linclass(data_file,W,[])
 *    [score,true_lab] = svmlight_linclass(data_file,W,W0)
 *    [score,true_lab] = svmlight_linclass(data_file,W,W0,verb)
 *
 * 
 *  Input:
 *   data_file [string] Path to file with examples stored in SVM^light format.
 *   W [nDims x nModels] Parameter vectors of nModels linear classifiers.
 *   W0 [nModels x 1] Bias of decision rule. If W0 is empty then W0 = 0 is used. 
 *   verb [1x1] If ~= 0 then prints info (default 0).
 *
 *  Output:
 *   score [nModels x nExamples] score(i,j) = W(:,i)'*X_j + W0(i)
 *   true_labels [nExamples x 1] labels from the data_file
 *
 *=================================================================*/

#define _FILE_OFFSET_BITS  64

#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <mex.h>
#include <sys/time.h>
#include <time.h>
#include <errno.h>

#include "libocas.h"
#include "ocas_helper.h"
#include "lib_svmlight_format.h"

#if !defined(MX_API_VER) || MX_API_VER<0x07040000
#define mwSize int
#define mwIndex int
#endif


/*----------------------------------------------------------
It computes dot produt between col-th column of matrix 
and a sparse vector given by vec_val and vec_idx (indices are 1 based) 
its length is vec_len. The matrix can be either Matlab full or Matlab 
sparse matrix.
-----------------------------------------------------------*/

double my_dot_prod(mxArray *matrix, uint32_t col, double *vec_val, uint32_t *vec_idx, long vec_len)
{

  double dot_prod;
  long i;


  if(!mxIsSparse(matrix)) {

    long mat_dims = mxGetM(matrix);
    double *mat_ptr = (double*)mxGetPr(matrix) + col*mat_dims;

    dot_prod = 0;
    for(i=0; i < vec_len; i++) 
      dot_prod += *(mat_ptr+vec_idx[i]-1) * vec_val[i];

  }
  else
  {

    uint32_t col_len, ptr1, ptr2, row1, row2;
    mwSize *Ir, *Jc;
    double *Pr, val1, val2;

    Ir = mxGetIr(matrix);
    Jc = mxGetJc(matrix);
    Pr = mxGetPr(matrix);

    col_len = Jc[col+1] - Jc[col];

    ptr1 = Jc[col];

    ptr2 = 0;

    dot_prod = 0;

    while(ptr2 < vec_len && (ptr1 - Jc[col]) < col_len) {
      row1 = Ir[ptr1];
      row2 = vec_idx[ptr2]-1;

      if(row1 == row2 ) {
        dot_prod += Pr[ptr1]*vec_val[ptr2];
        ptr1++;
        ptr2++;
      }
      else if(row1 < row2 ) {
        ptr1++;
      }
      else {
        ptr2 ++;
      }
    }

  }

  return(dot_prod);
}

/*======================================================================
  Main code plus interface to Matlab.
========================================================================*/

void mexFunction( int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[] )
{
  char fname[200];
  FILE *fid;
  char *line;
  double *feat_val;
  uint32_t *feat_idx;
  uint32_t nData = 0, nY = 0;
  long nnzf;
  uint64_t nnz = 0;
  mxArray *W;
  int nModels;
  double *score;
  double *W0;
  double *error;
  long max_dim = 0;
  long nDims = 0;
  long i, j;
  double *true_lab;
  double dot_prod;
  int verb = 0;
  double *auc;

  if( nrhs < 3 )
    mexErrMsgTxt("At least three input arguments must be passed.");

  if( nrhs > 4)
    mexErrMsgTxt("At most four input arguments must be passed.");

  if( nrhs == 4)
    verb = (int)mxGetScalar(prhs[3]);

  /* get input arguments */
  mxGetString(prhs[0], fname, 200);

  W = (mxArray*)prhs[1];
  nDims = mxGetM(prhs[1]);
  nModels = mxGetN(prhs[1]);

  if( mxIsEmpty(prhs[2]) == false)
  {
    W0 = (double*)mxGetPr(prhs[2]);
    if((LIBOCAS_MAX(mxGetM(prhs[2]),mxGetN(prhs[2])) != nModels) || 
       (LIBOCAS_MIN(mxGetM(prhs[2]),mxGetN(prhs[2])) != 1))
    {
      mexErrMsgTxt("The third input argument is of incorrect dimensionality.");
    }
  }
  else
  {
    W0 = mxCalloc(nModels, sizeof(double));
    if( W0 == NULL )
      mexErrMsgTxt("Not enough memmory to allocate W0.");
  }

  if(verb)
  {
    mexPrintf("input_file: %s\n", fname);
    mexPrintf("nModels: %d\n", nModels);
    mexPrintf("W dimensionality: %d\n", nDims);
  }

  fid = fopen(fname, "r");
  if(fid == NULL) {
    perror("fopen error: ");
    mexErrMsgTxt("Cannot open input file.");
  }

  /**********************************/
  line = mxCalloc(LIBSLF_MAXLINELEN, sizeof(char));
  if( line == NULL )
    mexErrMsgTxt("Not enogh memmory to allocate line buffer.");

  feat_idx = mxCalloc(LIBSLF_MAXLINELEN, sizeof(long));
  if( feat_idx == NULL )
    mexErrMsgTxt("Not enogh memmory to allocate feature buffer.");

  feat_val = mxCalloc(LIBSLF_MAXLINELEN, sizeof(double));
  if( feat_val == NULL )
    mexErrMsgTxt("Not enough memmory to allocate feature value buffer.");


  /*---------------------------------------------------------------------- */
  /* FIRST STEP: go through the data to get their number and dimensonality */
  /*---------------------------------------------------------------------- */

  if(verb) 
    mexPrintf("Analysing examples...");
  int label;
  int go = 1;
  long line_cnt = 0;

  while(go) {
    
    if(fgets(line,LIBSLF_MAXLINELEN, fid) == NULL ) 
    {
      go = 0;
      if(verb)
      {
        if( (line_cnt % 1000) != 0) 
          mexPrintf(" %ld", line_cnt);
        mexPrintf(" EOF.\n");
      }

    }
    else
    {
      line_cnt ++;
      nnzf = svmlight_format_parse_line(line, &label, feat_idx, feat_val);
      
      if(nnzf == -1) 
      {
           mexPrintf("Parsing error on line %ld .\n", line_cnt);
           mexErrMsgTxt("Probably defective input file.\n");
      }

      max_dim = LIBOCAS_MAX(max_dim,feat_idx[nnzf-1]);
      nnz += nnzf;
      nY = LIBOCAS_MAX(nY,label);
      
      if(verb && (line_cnt % 1000) == 0) {
        mexPrintf(" %ld", line_cnt);
        fflush(NULL);
      }
    }
  }

  /* rewind the file at the begining */
  fseek(fid, 0L, SEEK_SET);

  nData = line_cnt;
  if(verb)
  {
    mexPrintf("Data statistics:\n");
    mexPrintf("# of examples: %ld\n", nData);
    mexPrintf("# of classes : %ld\n", nY);
    mexPrintf("dimensionality: %d\n", max_dim);
    mexPrintf("nnz: %ld, density: %f%%\n", (long)nnz, 100*(double)nnz/((double)max_dim*(double)nData));
  }

  if(max_dim > nDims) 
  {
     mexErrMsgTxt("\nDimansionality of examples in the file exceeds dimensionality of W.");
  }    

  plhs[0] = (mxArray*)mxCreateDoubleMatrix(nModels,line_cnt,mxREAL);
  if(plhs[0] == NULL)
    mexErrMsgTxt("Not enought memory to allocate buffer for score.");
  score = (double*)mxGetPr(plhs[0]);

  plhs[1] = (mxArray*)mxCreateDoubleMatrix(line_cnt,1,mxREAL);
  if(plhs[1] == NULL)
    mexErrMsgTxt("Not enought memory to allocate buffer for true_labels.");
  true_lab = (double*)mxGetPr(plhs[1]);

  /*********************************************/
  /* MAIN LOOP                                 */
  /*********************************************/

  if(verb) mexPrintf("Reading examples...");

  go = 1;
  line_cnt = 0;

  while(go) {
    
    if(fgets(line,LIBSLF_MAXLINELEN, fid) == NULL ) 
    {
      go = 0;
      if(verb)
      {
        if( (line_cnt % 1000) != 0) 
          mexPrintf(" %ld", line_cnt);
        mexPrintf(" EOF.\n");
      }

    }
    else
    {
      line_cnt ++;
      nnzf = svmlight_format_parse_line(line, &label,feat_idx, feat_val);
      
      if(nnzf == -1) 
      {
         mexPrintf("Parsing error on line %ld .\n", line_cnt);
         mexErrMsgTxt("Probably defective input file.\n");
      }

      true_lab[line_cnt-1] = (double)label;

      for(i=0; i< nModels; i++) {

        dot_prod = my_dot_prod(W,i,feat_val,feat_idx,nnzf);
        score[LIBOCAS_INDEX(i,line_cnt-1,nModels)] = dot_prod + W0[i]; 
      } 

      if(verb && (line_cnt % 1000) == 0) {
        mexPrintf(" %ld", line_cnt);
        fflush(NULL);
      }

    }
  }

  fclose(fid);
  
  return;
}

