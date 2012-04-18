/*-----------------------------------------------------------------------
 * ocas_helper.c: Implementation of helper functions for the OCAS solver.
 *
 * It supports both sparse and dense matrices and loading data from
 * the SVM^light format.
 *-------------------------------------------------------------------- */

#define _FILE_OFFSET_BITS  64

#include <pthread.h>

#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <sys/time.h>
#include <stdlib.h>
#include <time.h>

#include "lib_svmlight_format.h"
#include "libocas.h"
#include "ocas_helper.h"

mxArray *data_X;
uint32_t nDim, nData, nY;
double *data_y;
cutting_plane_buf_T sparse_A;
double *full_A;
double *W;
double *oldW;
double *new_a;

double *A0;
double W0;
double oldW0;
double X0;

/* parallelization via threads */
struct thread_params_output
{
	double* output;
	uint32_t start;
	uint32_t end;
};

struct thread_qsort
{
	double* output;
/*	uint32_t* index;*/
	double* data;
	uint32_t size;
};

struct thread_params_add
{
  double *new_a;
  uint32_t *new_cut;
  uint32_t start;
  uint32_t end;
};


typedef enum 
{ 
  FALSE = 0,
  TRUE = 1
} 
boolean;


static int qsort_threads;
static pthread_t* threads = NULL;
static uint32_t* thread_slices = NULL;
static int num_threads;
static const int sort_limit=4096;
static struct thread_params_output* params_output;
static struct thread_params_add* params_add;

/* use multi-threads only if minimal number of examples to add is higher than the constant*/
static const uint32_t MinimalParallelCutLenght = 100;  

/***********************************************************************
   Multiclass SVM helper functions.                                   
 **********************************************************************/

/*----------------------------------------------------------------------
  sq_norm_W = sparse_compute_W( alpha, nSel ) does the following:

  oldW = W;
  W = sparse_A(:,1:nSel)'*alpha;
  sq_norm_W = W'*W;
  dp_WoldW = W'*oldW';

  ----------------------------------------------------------------------*/
void msvm_sparse_compute_W( double *sq_norm_W, 
                       double *dp_WoldW, 
                       double *alpha, 
                       uint32_t nSel, 
                       void* user_data )
{
  uint32_t i,j, nz_dims;

  memcpy(oldW, W, sizeof(double)*nY*nDim ); 
  memset(W, 0, sizeof(double)*nY*nDim );

  for(i=0; i < nSel; i++) {
    nz_dims = sparse_A.nz_dims[i];
    if(nz_dims > 0 && alpha[i] > 0) {
      for(j=0; j < nz_dims; j++) {
        W[sparse_A.index[i][j]] += alpha[i]*sparse_A.value[i][j];
      }
    }
  }

  *sq_norm_W = 0;
  *dp_WoldW = 0;
  for(j=0; j < nY*nDim; j++)
  {
    *sq_norm_W += W[j]*W[j];
    *dp_WoldW += W[j]*oldW[j];
  }
  
  return;
}


/*----------------------------------------------------------------------------------
  sq_norm_W = sparse_update_W( t ) does the following:

  W = oldW*(1-t) + t*W;
  sq_norm_W = W'*W;

  ---------------------------------------------------------------------------------*/
double msvm_full_update_W( double t, void* user_data )
{
  uint32_t j;
  double sq_norm_W;         

  sq_norm_W = 0;

  for(j=0; j < nY*nDim; j++) {
    W[j] = oldW[j]*(1-t) + t*W[j];
    sq_norm_W += W[j]*W[j];
  }          

  return( sq_norm_W );
}


/*----------------------------------------------------------------------------------
  sparse_add_new_cut( new_col_H, new_cut, cut_length, nSel ) does the following:

    new_a = zeros(nDim,nY);
    for i=1:nData
       if new_cut(i) ~= data_y(i)
          new_a(:,data_y(i)) = new_a(:,data_y(i)) + X(;,i);
          new_a(:,new_cut(i)) = new_a(:,new_cut(i)) - X(;,i);
       end
    end

    new_col_H = [sparse_A(:,1:nSel)'*new_a ; new_a'*new_a];
    sparse_A(:,nSel+1) = new_a;

    Warning: data_y is 1-based while new_cut is 0-based

  ---------------------------------------------------------------------------------*/
int msvm_sparse_add_new_cut( double *new_col_H, 
                         uint32_t *new_cut, 
                         uint32_t nSel, 
                         void* user_data )
{
/*  double *new_a, */
  double sq_norm_a;
  uint32_t i, j, nz_dims, ptr, y;

  memset(new_a, 0, sizeof(double)*nY*nDim);
  
  for(i=0; i < nData; i++)
  {
    y = (uint32_t)(data_y[i]-1);
    if(new_cut[i] != y)
    {
      add_sparse_col(&new_a[nDim*y], data_X, i);
      subtract_sparse_col(&new_a[nDim*(uint32_t)new_cut[i]], data_X, i);
    }
  }
 
  /* compute new_a'*new_a and count number of non-zero dimensions */
  nz_dims = 0; 
  sq_norm_a = 0;
  for(j=0; j < nY*nDim; j++ ) {
    if(new_a[j] != 0) {
      nz_dims++;
      sq_norm_a += new_a[j]*new_a[j];
    }
  }

  /* sparsify new_a and insert it to the last column  of sparse_A */
  sparse_A.nz_dims[nSel] = nz_dims;
  if(nz_dims > 0) {
    sparse_A.index[nSel] = NULL;
    sparse_A.value[nSel] = NULL;
    sparse_A.index[nSel] = mxCalloc(nz_dims,sizeof(uint32_t));
    sparse_A.value[nSel] = mxCalloc(nz_dims,sizeof(double));
    if(sparse_A.index[nSel]==NULL || sparse_A.value[nSel]==NULL)
    {
/*      mexErrMsgTxt("Not enough memory for vector sparse_A.index[nSel], sparse_A.value[nSel].");*/
      mxFree(sparse_A.index[nSel]);
      mxFree(sparse_A.value[nSel]);
      return(-1);
    }

    ptr = 0;
    for(j=0; j < nY*nDim; j++ ) {
      if(new_a[j] != 0) {
        sparse_A.index[nSel][ptr] = j;
        sparse_A.value[nSel][ptr++] = new_a[j];
      }
    }
  }
   
  new_col_H[nSel] = sq_norm_a;
  for(i=0; i < nSel; i++) {
    double tmp = 0;

    for(j=0; j < sparse_A.nz_dims[i]; j++) {
      tmp += new_a[sparse_A.index[i][j]]*sparse_A.value[i][j];
    }
      
    new_col_H[i] = tmp;
  }

  return 0;
}


/*----------------------------------------------------------------------
  sparse_compute_output( output ) does the follwing:

  output = W'*data_X;
  ----------------------------------------------------------------------*/
int msvm_sparse_compute_output( double *output, void* user_data )
{
  uint32_t i,y;

  for(i=0; i < nData; i++) 
  {
    for(y=0; y < nY; y++)
    {
      output[LIBOCAS_INDEX(y,i,nY)] = dp_sparse_col(&W[y*nDim], data_X, i);
    }
  }
  
  return 0;
}


/*-----------------------------------------------------------
  Functions working with full data.
  -------------------------------------------------------------*/

/*----------------------------------------------------------------------------------
  full_add_new_cut( new_col_H, new_cut, cut_length, nSel ) does the following:

    new_a = sum(data_X(:,find(new_cut ~=0 )),2);
    new_col_H = [full_A(:,1:nSel)'*new_a ; new_a'*new_a];
    full_A(:,nSel+1) = new_a;

  ---------------------------------------------------------------------------------*/
int msvm_full_add_new_cut( double *new_col_H, uint32_t *new_cut, uint32_t nSel, void* user_data)
{
  double sq_norm_a, *ptr;
  uint32_t i, j, y, y2;

  ptr = mxGetPr(data_X);

  memset(new_a, 0, sizeof(double)*nDim*nY);

  for(i=0; i < nData; i++)
  {
    y = (uint32_t)(data_y[i]-1);
    y2 = (uint32_t)new_cut[i];
    if(y2 != y)
    {
      for(j=0; j < nDim; j++ ) 
      {
        new_a[LIBOCAS_INDEX(j,y,nDim)] += ptr[LIBOCAS_INDEX(j,i,nDim)];
        new_a[LIBOCAS_INDEX(j,y2,nDim)] -= ptr[LIBOCAS_INDEX(j,i,nDim)];
      }
    }
  }

  /* compute new_a'*new_a and insert new_a to the last column of full_A */
  sq_norm_a = 0;
  for(j=0; j < nDim*nY; j++ ) {
    sq_norm_a += new_a[j]*new_a[j];
    full_A[LIBOCAS_INDEX(j,nSel,nDim*nY)] = new_a[j];
  }

  new_col_H[nSel] = sq_norm_a;
  for(i=0; i < nSel; i++) {
    double tmp = 0;

    for(j=0; j < nDim*nY; j++ ) {
      tmp += new_a[j]*full_A[LIBOCAS_INDEX(j,i,nDim*nY)];
    }
    new_col_H[i] = tmp;
  }

  return 0;
}


/*----------------------------------------------------------------------
  full_compute_output( output ) does the follwing:

  output = data_X'*W;
  ----------------------------------------------------------------------*/
int msvm_full_compute_output( double *output, void* user_data )
{
  uint32_t i, j, y;
  double *ptr, tmp;

  ptr = mxGetPr( data_X );

  for(i=0; i < nData; i++) 
  { 
    for(y=0; y < nY; y++)
    {
      tmp = 0;

      for(j=0; j < nDim; j++ ) 
      {
        tmp += W[LIBOCAS_INDEX(j,y,nDim)]*ptr[LIBOCAS_INDEX(j,i,nDim)];
      }
      
      output[LIBOCAS_INDEX(y,i,nY)] = tmp;
    }
  }
  
  return 0;
}



/*----------------------------------------------------------------------
  sq_norm_W = full_compute_W( alpha, nSel ) does the following:

  oldW = W;
  W = full_A(:,1:nSel)'*alpha;
  sq_norm_W = W'*W;
  dp_WoldW = W'*oldW';

  ----------------------------------------------------------------------*/
void msvm_full_compute_W( double *sq_norm_W, double *dp_WoldW, double *alpha, uint32_t nSel, void* user_data )
{
  uint32_t i,j;

  memcpy(oldW, W, sizeof(double)*nDim*nY ); 
  memset(W, 0, sizeof(double)*nDim*nY);

  for(i=0; i < nSel; i++) {
    if( alpha[i] > 0 ) {
      for(j=0; j< nDim*nY; j++ ) {
        W[j] += alpha[i]*full_A[LIBOCAS_INDEX(j,i,nDim*nY)];
      }

    }
  }

  *sq_norm_W = 0;
  *dp_WoldW = 0;
  for(j=0; j < nDim*nY; j++) {
    *sq_norm_W += W[j]*W[j];
    *dp_WoldW += W[j]*oldW[j];
  }

  return;
}



/***********************************************************************
  Generic functions.
 **********************************************************************/


/*-----------------------------------------------------------------------
  Print statistics.
  -----------------------------------------------------------------------*/
void ocas_print(ocas_return_value_T value)
{
  mexPrintf("%4d: tim=%f, Q_P=%f, Q_D=%f, Q_P-Q_D=%f, 1-Q_D/Q_P=%f, nza=%4d, err=%.2f%%, qpf=%d\n",
            value.nIter,value.ocas_time, value.Q_P,value.Q_D,value.Q_P-value.Q_D,(value.Q_P-value.Q_D)/LIBOCAS_ABS(value.Q_P), 
            value.nNZAlpha, 100*(double)value.trn_err/(double)nData, value.qp_exitflag );
}

void ocas_print_null(ocas_return_value_T value)
{
  return;
}


/*-----------------------------------------------------------------------
  Get absolute time in seconds.
  -----------------------------------------------------------------------*/
double get_time()
{
	struct timeval tv;
	if (gettimeofday(&tv, NULL)==0)
		return tv.tv_sec+((double)(tv.tv_usec))/1e6;
	else
		return 0.0;
}


/*=========================================================================
 *
 * Ocas helper functions implemented for input data represented as 
 * Matlab sparse matrix.
 *
 *=========================================================================*/

/*----------------------------------------------------------------------
  in-place computes sparse_mat(:,col)= alpha * sparse_mat(:,col)
  where alpha is a scalar and sparse_mat is Matlab sparse matrix.
  ----------------------------------------------------------------------*/
void mul_sparse_col(double alpha, mxArray *sparse_mat, uint32_t col)
{
	uint32_t nItems, ptr, i;
	INDEX_TYPE_T *Ir, *Jc;
	double *Pr;

	Ir = mxGetIr(sparse_mat);
	Jc = mxGetJc(sparse_mat);
	Pr = mxGetPr(sparse_mat);

	nItems = Jc[col+1] - Jc[col];
	ptr = Jc[col];

	for(i=0; i < nItems; i++)
		Pr[ptr++]*=alpha;
}


/*----------------------------------------------------------------------
 It computes full_vec = full_vec + sparse_mat(:,col)
 where full_vec is a double array and sparse_mat is Matlab 
 sparse matrix.
  ----------------------------------------------------------------------*/
void add_sparse_col(double *full_vec, mxArray *sparse_mat, uint32_t col)
{
  uint32_t nItems, ptr, i, row;
  INDEX_TYPE_T *Ir, *Jc;
  double *Pr, val;
    
  Ir = mxGetIr(sparse_mat);
  Jc = mxGetJc(sparse_mat);
  Pr = mxGetPr(sparse_mat);

  nItems = Jc[col+1] - Jc[col];
  ptr = Jc[col];

  for(i=0; i < nItems; i++) {
    val = Pr[ptr];
    row = Ir[ptr++];

    full_vec[row] += val;
  }
}

/*----------------------------------------------------------------------
 It computes full_vec = full_vec - sparse_mat(:,col)
 where full_vec is a double array and sparse_mat is Matlab 
 sparse matrix.
  ----------------------------------------------------------------------*/
void subtract_sparse_col(double *full_vec, mxArray *sparse_mat, uint32_t col)
{
  uint32_t nItems, ptr, i, row;
  INDEX_TYPE_T *Ir, *Jc;
  double *Pr, val;
    
  Ir = mxGetIr(sparse_mat);
  Jc = mxGetJc(sparse_mat);
  Pr = mxGetPr(sparse_mat);

  nItems = Jc[col+1] - Jc[col];
  ptr = Jc[col];

  for(i=0; i < nItems; i++) {
    val = Pr[ptr];
    row = Ir[ptr++];

    full_vec[row] -= val;
  }
}

/*----------------------------------------------------------------------
 It computes dp = full_vec'*sparse_mat(:,col)
 where full_vec is a double array and sparse_mat is Matlab 
 sparse matrix.
  ----------------------------------------------------------------------*/
double dp_sparse_col(double *full_vec, mxArray *sparse_mat, uint32_t col)
{
  uint32_t nItems, ptr, i, row;
  INDEX_TYPE_T *Ir, *Jc;
  double *Pr, val, dp;

  Ir = mxGetIr(sparse_mat);
  Jc = mxGetJc(sparse_mat);
  Pr = mxGetPr(sparse_mat);

  dp = 0;
  nItems = Jc[col+1] - Jc[col];
  ptr = Jc[col];

  for(i=0; i < nItems; i++) {
    val = Pr[ptr];
    row = Ir[ptr++];

    dp += full_vec[row]*val;
  }

  return(dp);  
}


/*----------------------------------------------------------------------------------
  sq_norm_W = sparse_update_W( t ) does the following:

  W = oldW*(1-t) + t*W;
  sq_norm_W = W'*W;

  ---------------------------------------------------------------------------------*/
double sparse_update_W( double t, void* user_data )
{
  uint32_t j;
  double sq_norm_W;         

  W0 = oldW0*(1-t) + t*W0;
  sq_norm_W = W0*W0;

  for(j=0; j <nDim; j++) {
    W[j] = oldW[j]*(1-t) + t*W[j];
    sq_norm_W += W[j]*W[j];
  }          

  return( sq_norm_W );
}

/*-------------------------------------------------------------------------
  sq_norm_W = full_update_W( t ) does the following:

  W = oldW*(1-t) + t*W;
  sq_norm_W = W'*W;
---------------------------------------------------------------------------*/
double full_update_W( double t, void* user_data )
{
  uint32_t j;
  double sq_norm_W;         

  W0 = oldW0*(1-t) + t*W0;
  sq_norm_W = W0*W0;

  for(j=0; j <nDim; j++) {
    W[j] = oldW[j]*(1-t) + t*W[j];
    sq_norm_W += W[j]*W[j];
  }          

  return( sq_norm_W );
}


/*----------------------------------------------------------------------------------
  sparse_add_new_cut( new_col_H, new_cut, cut_length, nSel ) does the following:

    new_a = sum(data_X(:,find(new_cut ~=0 )),2);
    new_col_H = [sparse_A(:,1:nSel)'*new_a ; new_a'*new_a];
    sparse_A(:,nSel+1) = new_a;

  ---------------------------------------------------------------------------------*/
int sparse_add_new_cut( double *new_col_H, 
                         uint32_t *new_cut, 
                         uint32_t cut_length, 
                         uint32_t nSel, 
                         void* user_data )
{
/*  double *new_a, */
  double sq_norm_a;
  uint32_t i, j, nz_dims, ptr;

  memset(new_a, 0, sizeof(double)*nDim);
  
  for(i=0; i < cut_length; i++) {
    add_sparse_col(new_a, data_X, new_cut[i]);

    A0[nSel] += X0*data_y[new_cut[i]];    
  }
 
  /* compute new_a'*new_a and count number of non-zero dimensions */
  nz_dims = 0; 
  sq_norm_a = A0[nSel]*A0[nSel];
  for(j=0; j < nDim; j++ ) {
    if(new_a[j] != 0) {
      nz_dims++;
      sq_norm_a += new_a[j]*new_a[j];
    }
  }

  /* sparsify new_a and insert it to the last column  of sparse_A */
  sparse_A.nz_dims[nSel] = nz_dims;
  if(nz_dims > 0) {
    sparse_A.index[nSel] = NULL;
    sparse_A.value[nSel] = NULL;
    sparse_A.index[nSel] = mxCalloc(nz_dims,sizeof(uint32_t));
    sparse_A.value[nSel] = mxCalloc(nz_dims,sizeof(double));
    if(sparse_A.index[nSel]==NULL || sparse_A.value[nSel]==NULL)
    {
/*      mexErrMsgTxt("Not enough memory for vector sparse_A.index[nSel], sparse_A.value[nSel].");*/
      mxFree(sparse_A.index[nSel]);
      mxFree(sparse_A.value[nSel]);
      return(-1);
    }

    ptr = 0;
    for(j=0; j < nDim; j++ ) {
      if(new_a[j] != 0) {
        sparse_A.index[nSel][ptr] = j;
        sparse_A.value[nSel][ptr++] = new_a[j];
      }
    }
  }
   
  new_col_H[nSel] = sq_norm_a;
  for(i=0; i < nSel; i++) {
    double tmp = A0[nSel]*A0[i];

    for(j=0; j < sparse_A.nz_dims[i]; j++) {
      tmp += new_a[sparse_A.index[i][j]]*sparse_A.value[i][j];
    }
      
    new_col_H[i] = tmp;
  }

/*  mxFree( new_a );*/

  return 0;
}


/*----------------------------------------------------------------------------------
  full_add_new_cut( new_col_H, new_cut, cut_length, nSel ) does the following:

    new_a = sum(data_X(:,find(new_cut ~=0 )),2);
    new_col_H = [full_A(:,1:nSel)'*new_a ; new_a'*new_a];
    full_A(:,nSel+1) = new_a;

  ---------------------------------------------------------------------------------*/
int full_add_new_cut( double *new_col_H, 
                       uint32_t *new_cut, 
                       uint32_t cut_length, 
                       uint32_t nSel,
                       void* user_data)
{
/*  double *new_a, */
  double sq_norm_a, *ptr;
  uint32_t i, j;

  ptr = mxGetPr(data_X);

  memset(new_a, 0, sizeof(double)*nDim);


  for(i=0; i < cut_length; i++) {
    for(j=0; j < nDim; j++ ) {
      new_a[j] += ptr[LIBOCAS_INDEX(j,new_cut[i],nDim)];
    }

    A0[nSel] += X0*data_y[new_cut[i]];    
  }

  /* compute new_a'*new_a and insert new_a to the last column of full_A */
  sq_norm_a = A0[nSel]*A0[nSel];
  for(j=0; j < nDim; j++ ) {
    sq_norm_a += new_a[j]*new_a[j];
    full_A[LIBOCAS_INDEX(j,nSel,nDim)] = new_a[j];
  }

  new_col_H[nSel] = sq_norm_a;
  for(i=0; i < nSel; i++) {
    double tmp = A0[nSel]*A0[i];

    for(j=0; j < nDim; j++ ) {
      tmp += new_a[j]*full_A[LIBOCAS_INDEX(j,i,nDim)];
    }
    new_col_H[i] = tmp;
  }

/*  mxFree( new_a );*/

  return 0;
}


/*----------------------------------------------------------------------
  sparse_compute_output( output ) does the follwing:

  output = data_X'*W;
  ----------------------------------------------------------------------*/
int sparse_compute_output( double *output, void* user_data )
{
  uint32_t i;

  for(i=0; i < nData; i++) { 
    output[i] = data_y[i]*X0*W0 + dp_sparse_col(W, data_X, i);
  }
  
  return 0;
}

/*----------------------------------------------------------------------
  full_compute_output( output ) does the follwing:

  output = data_X'*W;
  ----------------------------------------------------------------------*/
int full_compute_output( double *output, void* user_data )
{
  uint32_t i, j;
  double *ptr, tmp;

  ptr = mxGetPr( data_X );

  for(i=0; i < nData; i++) { 
    tmp = data_y[i]*X0*W0;

    for(j=0; j < nDim; j++ ) {
      tmp += W[j]*ptr[LIBOCAS_INDEX(j,i,nDim)];
    }
    output[i] = tmp;
  }
  
  return 0;
}



/*----------------------------------------------------------------------
  sq_norm_W = sparse_compute_W( alpha, nSel ) does the following:

  oldW = W;
  W = sparse_A(:,1:nSel)'*alpha;
  sq_norm_W = W'*W;
  dp_WoldW = W'*oldW';

  ----------------------------------------------------------------------*/
void sparse_compute_W( double *sq_norm_W, 
                       double *dp_WoldW, 
                       double *alpha, 
                       uint32_t nSel, 
                       void* user_data )
{
  uint32_t i,j, nz_dims;

  memcpy(oldW, W, sizeof(double)*nDim ); 
  memset(W, 0, sizeof(double)*nDim);

  oldW0 = W0;
  W0 = 0;

  for(i=0; i < nSel; i++) {
    nz_dims = sparse_A.nz_dims[i];
    if(nz_dims > 0 && alpha[i] > 0) {
      for(j=0; j < nz_dims; j++) {
        W[sparse_A.index[i][j]] += alpha[i]*sparse_A.value[i][j];
      }
    }
    W0 += A0[i]*alpha[i];
  }

  *sq_norm_W = W0*W0;
  *dp_WoldW = W0*oldW0;
  for(j=0; j < nDim; j++) {
    *sq_norm_W += W[j]*W[j];
    *dp_WoldW += W[j]*oldW[j];
  }
  
  return;
}


/*----------------------------------------------------------------------
  sq_norm_W = full_compute_W( alpha, nSel ) does the following:

  oldW = W;
  W = full_A(:,1:nSel)'*alpha;
  sq_norm_W = W'*W;
  dp_WoldW = W'*oldW';

  ----------------------------------------------------------------------*/
void full_compute_W( double *sq_norm_W, double *dp_WoldW, double *alpha, uint32_t nSel, void* user_data )
{
  uint32_t i,j;

  memcpy(oldW, W, sizeof(double)*nDim ); 
  memset(W, 0, sizeof(double)*nDim);

  oldW0 = W0;
  W0 = 0;

  for(i=0; i < nSel; i++) {
    if( alpha[i] > 0 ) {
      for(j=0; j< nDim; j++ ) {
        W[j] += alpha[i]*full_A[LIBOCAS_INDEX(j,i,nDim)];
      }

      W0 += A0[i]*alpha[i];
    }
  }

  *sq_norm_W = W0*W0;
  *dp_WoldW = W0*oldW0;
  for(j=0; j < nDim; j++) {
    *sq_norm_W += W[j]*W[j];
    *dp_WoldW += W[j]*oldW[j];
  }

  return;
}


/* ==========================================================================
 Multi-thread version of the OCAS helper functions.

 int parallel_sparse_compute_output( double *output, void* user_data )

 void destroy_parallel_ocas(void)
 int init_parallel_ocas(int number_of_threads)

============================================================================*/


/*----------------------------------------------------------------------
  parallel_sparse_compute_output( output, user_data ) does the following

  output = data_X'*W + W0*X0*y[i];
  ----------------------------------------------------------------------*/

static void* parallel_compute_output_helper(void* p)
{
	struct thread_params_output* params = (struct thread_params_output*) p;
	double* output=params->output;
	uint32_t start=params->start;
	uint32_t end=params->end;
	uint32_t i;

    for(i=start; i < end; i++) 
      output[i] = data_y[i]*X0*W0 + dp_sparse_col(W, data_X, i);

    return(NULL);
}

int parallel_sparse_compute_output( double *output, void* user_data )
{
  /*  one-thraed code looks like:

  uint32_t i;

  for(i=0; i < nData; i++) { 
    output[i] = data_y[i]*X0*W0 + dp_sparse_col(W, data_X, i);
  }
  */

/*  struct thread_params_output params;*/

  int nthreads=num_threads-1;
  int end=0;
  int t;
  
  if (nData < num_threads)
  {
    nthreads=nData-1;
  }

  for (t=0; t<nthreads; t++)
  {
    params_output[t].output = output;
    if (t==0)
      params_output[t].start = 0;
    else
      params_output[t].start = thread_slices[t-1];
    params_output[t].end = thread_slices[t];
    
    if (pthread_create(&threads[t], NULL, parallel_compute_output_helper, (void*)&params_output[t]) != 0)
    {
      mexPrintf("\nError: Thread creation failed.\n");
      return(-1);
    }

    end=params_output[t].end;
  }

  params_output[t].output = output;
  params_output[t].start = end;
  params_output[t].end = nData;

  parallel_compute_output_helper(&params_output[t]);

  for (t=0; t<nthreads; t++)
  {
    if (pthread_join(threads[t], NULL) != 0)
    {
      mexPrintf("\nError: pthread_join failed.\n");
      return(-1);
    }
  }

  return 0;
}

/* Ihis function initialize parallel processing . It must be
 called prior to paralel_XXXX  OCAS heleper functions. */
int init_parallel_ocas(int number_of_threads)
{
  num_threads = number_of_threads;

  thread_slices = (uint32_t*)mxCalloc(num_threads,sizeof(uint32_t));
  if(thread_slices == NULL) 
  {
    mexPrintf("Not enough memory for vector num_threads.");
    goto clean_up;

  }

  threads = (pthread_t*)mxCalloc(num_threads,sizeof(pthread_t));
  if(threads== NULL) 
  {
    mexPrintf("Not enough memory for threads structure.");
    goto clean_up;
  }

  params_output = (struct thread_params_output*)mxCalloc(num_threads,sizeof(struct thread_params_output));
  if(params_output== NULL) 
  {
    mexPrintf("Not enough memory for params structure.");
    goto clean_up;
  }

  params_add = (struct thread_params_add*)mxCalloc(num_threads,sizeof(struct thread_params_add));
  if(params_add== NULL) 
  {
    mexPrintf("Not enough memory for params structure.");
    goto clean_up;
  }

  /* The following code finds splits of the data such that each
  data slice contains approximately the same number of nonzero elements. */ 
  uint32_t i;
  int nnz_split = 0;
  INDEX_TYPE_T* Jc = mxGetJc(data_X);

  for(i=0; i < nData; i++)
    nnz_split += Jc[i+1] - Jc[i];

  nnz_split/=num_threads;
  uint64_t accum_nnz = 0;
  int thr = 0;

  mexPrintf("Data slices: ");
  for(i=0; i < nData; i++)
  {
    INDEX_TYPE_T nItems = Jc[i+1] - Jc[i];

    if (accum_nnz < nnz_split*(thr+1))
      accum_nnz+=nItems;
    else
    {
      thread_slices[thr]=i;
      mexPrintf("%d ", thread_slices[thr]);
      accum_nnz+=nItems;
      thr++;
    }
  }
  mexPrintf("\n");
  
  return(0);

clean_up:    

  mxFree(threads);
  mxFree(thread_slices);

  return(-1);
}

/* release memory allocated for parallelized functions */
void destroy_parallel_ocas(void)
{
  mxFree(params_add);
  mxFree(params_output);
  mxFree(threads);
  mxFree(thread_slices);
}

/*----------------------------------------------------------------------------------
  Parallel version of 
  sparse_add_new_cut( new_col_H, new_cut, cut_length, nSel ) does the following:

    new_a = sum(data_X(:,find(new_cut ~=0 )),2);
    new_col_H = [sparse_A(:,1:nSel)'*new_a ; new_a'*new_a];
    sparse_A(:,nSel+1) = new_a;

  ---------------------------------------------------------------------------------*/

static void* parallel_sparse_add_helper(void* p)
{
	struct thread_params_add* params = (struct thread_params_add*) p;
	double* local_new_a=params->new_a;
    uint32_t* new_cut=params->new_cut;
	uint32_t start=params->start;
	uint32_t end=params->end;
	uint32_t i;

    for(i=start; i <= end; i++) 
      add_sparse_col(local_new_a, data_X, new_cut[i]);

    return(NULL);
}

int parallel_sparse_add_new_cut( double *new_col_H, 
                         uint32_t *new_cut, 
                         uint32_t cut_length, 
                         uint32_t nSel, 
                         void* user_data )
{
/*  double *new_a, */
  double sq_norm_a;
  uint32_t i, j, nz_dims, ptr;

  /* temporary vector */
/*  new_a = (double*)mxCalloc(nDim,sizeof(double));*/
/*  if(new_a == NULL) */
/*    return(-1); */
  memset(new_a, 0, sizeof(double)*nDim*num_threads);
  
  if((cut_length < MinimalParallelCutLenght) && (cut_length >= num_threads))
  {

    for(i=0; i < cut_length; i++) {
      add_sparse_col(new_a, data_X, new_cut[i]);

      A0[nSel] += X0*data_y[new_cut[i]];    
    }
  }
  else
  {
    
/*    struct thread_params_add params;*/
/*    params.new_cut = new_cut;*/

    uint32_t chunk = cut_length/num_threads;
    uint32_t start = 0;
    uint32_t end = chunk-1;
    int t;
    for(t = 0; t < num_threads-1; t++)
    {
      params_add[t].start = start;
      params_add[t].end = end;
      params_add[t].new_a = &new_a[nDim*(t+1)];
      params_add[t].new_cut = new_cut;

      start = end+1;
      end = end+chunk;

      if (pthread_create(&threads[t], NULL, parallel_sparse_add_helper, (void*)&params_add[t]) != 0)
      {
        mexPrintf("\nError: Thread creation failed.\n");
        return(-1);
      }
       
    }

    params_add[t].start = start;
    params_add[t].end = cut_length-1;
    params_add[t].new_a = new_a;
    params_add[t].new_cut = new_cut;

    parallel_sparse_add_helper((void*)&params_add[t]);

    for (t=0; t<num_threads-1; t++)
    {
      if (pthread_join(threads[t], NULL) != 0)
      {
        return(-1);
      }

      double* a = &new_a[nDim*(t+1)];

      for (i=0; i<nDim; i++)
        new_a[i]+=a[i];
    }

    for(i=0; i < cut_length; i++) 
      A0[nSel] += X0*data_y[new_cut[i]];    
  }
 
  /* compute new_a'*new_a and count number of non-zero dimensions */
  nz_dims = 0; 
  sq_norm_a = A0[nSel]*A0[nSel];
  for(j=0; j < nDim; j++ ) {
    if(new_a[j] != 0) {
      nz_dims++;
      sq_norm_a += new_a[j]*new_a[j];
    }
  }

  /* sparsify new_a and insert it to the last column  of sparse_A */
  sparse_A.nz_dims[nSel] = nz_dims;
  if(nz_dims > 0) {
    sparse_A.index[nSel] = NULL;
    sparse_A.value[nSel] = NULL;
    sparse_A.index[nSel] = mxCalloc(nz_dims,sizeof(uint32_t));
    sparse_A.value[nSel] = mxCalloc(nz_dims,sizeof(double));
    if(sparse_A.index[nSel]==NULL || sparse_A.value[nSel]==NULL)
    {
/*      mexErrMsgTxt("Not enough memory for vector sparse_A.index[nSel], sparse_A.value[nSel].");*/
      mxFree(sparse_A.index[nSel]);
      mxFree(sparse_A.value[nSel]);
      return(-1);
    }

    ptr = 0;
    for(j=0; j < nDim; j++ ) {
      if(new_a[j] != 0) {
        sparse_A.index[nSel][ptr] = j;
        sparse_A.value[nSel][ptr++] = new_a[j];
      }
    }
  }
   
  new_col_H[nSel] = sq_norm_a;
  for(i=0; i < nSel; i++) {
    double tmp = A0[nSel]*A0[i];

    for(j=0; j < sparse_A.nz_dims[i]; j++) {
      tmp += new_a[sparse_A.index[i][j]]*sparse_A.value[i][j];
    }
      
    new_col_H[i] = tmp;
  }

/*  mxFree( new_a );*/

  return 0;
}


/*=========================================================================
 *
 * Ocas helper functions implemented for input data represented as
 * an array of doubles.
 *
 *=========================================================================*/



/*=======================================================================
 OCAS helper functions for sorting numbers.
=======================================================================*/
static void swapf(double* a, double* b)
{
	double dummy=*b;
	*b=*a;
	*a=dummy;
}

/*static void swapi(uint32_t* a, uint32_t* b)*/
/*{*/
/*	int dummy=*b;*/
/*	*b=*a;*/
/*	*a=dummy;*/
/*}*/

/* sort arrays value and data according to value in ascending order */
int qsort_data(double* value, double* data, uint32_t size)
{
    if(size == 1)
      return 0;

	if (size==2)
	{
		if (value[0] > value[1])
		{
			swapf(&value[0], &value[1]);
/*			swapi(&data[0], &data[1]);*/
			swapf(&data[0], &data[1]);
		}
		return 0;
	}
	double split=value[size/2];

	uint32_t left=0;
	uint32_t right=size-1;

	while (left<=right)
	{
		while (value[left] < split)
			left++;
		while (value[right] > split)
			right--;

		if (left<=right)
		{
			swapf(&value[left], &value[right]);
/*			swapi(&data[left], &data[right]);*/
			swapf(&data[left], &data[right]);
			left++;
			right--;
		}
	}

	if (right+1> 1)
		qsort_data(value,data,right+1);

	if (size-left> 1)
		qsort_data(&value[left],&data[left], size-left);


    return 0;
}

/*-------------------------------------------------------------------------
 Parallel version of qsort_data.
  -------------------------------------------------------------------------*/
void* parallel_qsort_helper(void* p)
{
	struct thread_qsort* ps=(struct thread_qsort*) p;
	double* output=ps->output;
/*	uint32_t* data=ps->data;*/
	double* data=ps->data;
	uint32_t size=ps->size;

    if(size == 1)
      return 0;


	if (size==2)
	{
		if (output[0] > output [1])
		{
			swapf(&output[0], &output[1]);
/*			swapi(&data[0], &data[1]);*/
			swapf(&data[0], &data[1]);
		}
		return(NULL);
	}
	/*double split=output[(((uint64_t) size)*rand())/(((uint64_t)RAND_MAX)+1)];*/
	double split=output[size/2];

	uint32_t left=0;
	uint32_t right=size-1;

	while (left<=right)
	{
		while (output[left] < split)
			left++;
		while (output[right] > split)
			right--;

		if (left<=right)
		{
			swapf(&output[left], &output[right]);
/*			swapi(&index[left], &index[right]);*/
			swapf(&data[left], &data[right]);
			left++;
			right--;
		}
	}
	boolean lthread_start=FALSE;
	boolean rthread_start=FALSE;
	pthread_t lthread;
	pthread_t rthread;
	struct thread_qsort t1;
	struct thread_qsort t2;

	if (right+1> 1 && (right+1< sort_limit || qsort_threads >= num_threads-1))
		qsort_data(output,data,right+1);
	else if (right+1> 1)
	{
		qsort_threads++;
		lthread_start=TRUE;
		t1.output=output;
		t1.data=data;
		t1.size=right+1;
		if (pthread_create(&lthread, NULL, parallel_qsort_helper, (void*) &t1) != 0)
		{
			lthread_start=FALSE;
			qsort_threads--;
			qsort_data(output,data,right+1);
		}
	}


	if (size-left> 1 && (size-left< sort_limit || qsort_threads >= num_threads-1))
		qsort_data(&output[left],&data[left], size-left);
	else if (size-left> 1)
	{
		qsort_threads++;
		rthread_start=TRUE;
		t2.output=&output[left];
		t2.data=&data[left];
		t2.size=size-left;
		if (pthread_create(&rthread, NULL, parallel_qsort_helper, (void*)&t2) != 0)
		{
			rthread_start=FALSE;
			qsort_threads--;
			qsort_data(&output[left],&data[left], size-left);
		}
	}

	if (lthread_start)
	{
		pthread_join(lthread, NULL);
		qsort_threads--;
	}

	if (rthread_start)
	{
		pthread_join(rthread, NULL);
		qsort_threads--;
	}

    return(NULL);
}


int parallel_qsort_data(double* value, double* data, uint32_t size)
{
  struct thread_qsort qthr;

  qsort_threads=0;

  qthr.output=value;
  qthr.data=data;
  qthr.size=size;
  parallel_qsort_helper((void*)&qthr);

  return(0);
}

/* =================================================================================
 Other auxiliary functions.
   =================================================================================*/

/* ---------------------------------------------------------------------------------
This function loads regularization constants from a text file. Each line contains
a single constant. 
  ---------------------------------------------------------------------------------*/
int load_regconsts(char *fname, double **vec_C, uint32_t *len_vec_C, int verb)
{
  double C;
  char *line = NULL;
  int exitflag = 0;
  FILE *fid;

  if(verb) mexPrintf("Input file: %s\n", fname);

  fid = fopen(fname, "r");
  if(fid == NULL) {
    perror("fopen error ");
    mexPrintf("Cannot open input file.\n");
    exitflag = -1;
    goto clean_up;
  }

  line = mxCalloc(LIBSLF_MAXLINELEN, sizeof(char));
  if( line == NULL )
  {
    mexPrintf("Not enough memmory to allocate line buffer.\n");
    exitflag = -1;
    goto clean_up;
  }


  if(verb) mexPrintf("Counting regularization constants...");
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

      C = atof(line);
      
      if(verb)
      {
        if( (line_cnt % 1000) == 0) {
          mexPrintf(" %ld", line_cnt);
          fflush(NULL);
        }
      }
    }
  }
  
  *vec_C = (double*)mxCalloc(line_cnt, sizeof(double));
  if( vec_C == NULL )
  {
    mexPrintf("Not enough memmory to allocate vec_C.\n");
    exitflag = -1;
    goto clean_up;
  }

  fseek(fid, 0, SEEK_SET);

  if(verb) mexPrintf("Reading regularization constants...");
  go = 1;
  line_cnt = 0;
  while(go) 
  {
    
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
      (*vec_C)[line_cnt] = atof(line);
      line_cnt ++;
      
      if(verb)
      {
        if( (line_cnt % 1000) == 0) {
          mexPrintf(" %ld", line_cnt);
          fflush(NULL);
        }
      }
    }
  }
  
  fclose(fid);
  *len_vec_C = line_cnt;

clean_up:
  mxFree(line);

  return(exitflag); 
}



/* --------------------------------------------------------------------------------- 
This function loads SVMlight data file to sparse matlab matrix data_X and 
dense vector data_y which both are assumed to global variables.
  ---------------------------------------------------------------------------------*/
int load_svmlight_file(char *fname, int verb)
{
  char *line = NULL;
  FILE *fid;
  double *feat_val = NULL;
  double sparse_memory_requirement, full_memory_requirement;
  uint32_t *feat_idx = NULL;
  long nnzf;
  int max_dim = 0;
  long j;
  uint64_t nnz = 0;
  mwSize *irs = NULL, *jcs = NULL;
  int exitflag = 0;
  double *sr = NULL;
  
/*  mexPrintf("Input file: %s\n", fname);*/

  fid = fopen(fname, "r");
  if(fid == NULL) {
    perror("fopen error ");
    mexPrintf("Cannot open input file.\n");
    exitflag = -1;
    goto clean_up;
  }

  line = mxCalloc(LIBSLF_MAXLINELEN, sizeof(char));
  if( line == NULL )
  {
    mexPrintf("Not enough memmory to allocate line buffer.\n");
    exitflag = -1;
    goto clean_up;
  }

  feat_idx = mxCalloc(LIBSLF_MAXLINELEN, sizeof(uint32_t));
  if( feat_idx == NULL )
  {
    mexPrintf("Not enough memmory to allocate feat_idx.\n");
    exitflag = -1;
    goto clean_up;
  }

  feat_val = mxCalloc(LIBSLF_MAXLINELEN, sizeof(double));
  if( feat_val == NULL )
  {
    mexPrintf("Not enough memmory to allocate feat_val.\n");
    exitflag = -1;
    goto clean_up;
  }

  if(verb) mexPrintf("Analysing input data...");
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
         mexPrintf("Probably defective input file.\n");
         exitflag = -1;
         goto clean_up;
      }

      max_dim = LIBOCAS_MAX(max_dim,feat_idx[nnzf-1]);
      nnz += nnzf;
      
      if(verb)
      {
        if( (line_cnt % 1000) == 0) {
          mexPrintf(" %ld", line_cnt);
          fflush(NULL);
        }
      }
    }
  }

  fclose(fid);  
  if(verb)
  {
    mexPrintf("Data statistics:\n");
    mexPrintf("# of examples: %ld\n", line_cnt);
    mexPrintf("dimensionality: %d\n", max_dim);
    mexPrintf("nnz: %ld, density: %f%%\n", (long)nnz, 100*(double)nnz/((double)max_dim*(double)line_cnt));
  }

  sparse_memory_requirement = ((double)nnz*((double)sizeof(double)+(double)sizeof(mwSize)))/(1024.0*1024.0);
  full_memory_requirement = sizeof(double)*(double)max_dim*(double)line_cnt/(1024.0*1024.0);

  if(verb)
  {
    mexPrintf("Memory requirements for sparse matrix: %.3f MB\n", sparse_memory_requirement);
    mexPrintf("Memory requirements for full matrix: %.3f MB\n", full_memory_requirement);
  }

  if( full_memory_requirement < sparse_memory_requirement)
  {
    if(verb)
      mexPrintf("Full matrix represenation used.\n");

    data_X = mxCreateDoubleMatrix(max_dim, line_cnt, mxREAL);

    if( data_X == NULL)
    {
      mexPrintf("Not enough memory to allocate data_X .\n");
      exitflag = -1;
      goto clean_up;
    }

  }
  else
  {
    if(verb)
      mexPrintf("Sparse matrix represenation used.\n");

    data_X = mxCreateSparse(max_dim, line_cnt, nnz, mxREAL);
    if( data_X == NULL)
    {
      mexPrintf("Not enough memory to allocate data_X .\n");
      exitflag = -1;
      goto clean_up;
    }

    sr  = mxGetPr(data_X);
    irs = (mwSize*)mxGetIr(data_X);
    jcs = (mwSize*)mxGetJc(data_X);

  }


/*  mexPrintf("Required memory: %.3f MB\n", */
/*    ((double)nnz*((double)sizeof(double)+(double)sizeof(mwSize)))/(1024.0*1024.0));*/

  /*---------------------------------------------*/

  data_y = mxCalloc(line_cnt, sizeof(double));
  if(data_y == NULL)
  {
    mexPrintf("Not enough memory to allocate data_y.\n");
    exitflag = -1;
    goto clean_up;
  }

  fid = fopen(fname, "r");
  if(fid == NULL) {
    perror("fopen error ");
    mexPrintf("Cannot open input file.\n");
    exitflag = -1;
    goto clean_up;
  }

  if(verb)
    mexPrintf("Reading examples...");
  
  go = 1;
  line_cnt = 0;
  long k=0;
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
         mexPrintf("Defective input file.\n");
         exitflag = -1;
         goto clean_up;
      }

      data_y[line_cnt-1] = (double)label;

      if( mxIsSparse( data_X) )
      {
        jcs[line_cnt-1] = k;

        for(j = 0; j < nnzf; j++) {
          sr[k] = feat_val[j];
          irs[k] = feat_idx[j]-1;
          k++;
        }
      }
      else
      {
        double *ptr = mxGetPr(data_X);
        for(j=0; j < nnzf; j++ ) {
          ptr[LIBOCAS_INDEX(feat_idx[j]-1,line_cnt-1,max_dim)] = feat_val[j];
        }

      }
      
      if(verb)
      {
        if( (line_cnt % 1000) == 0) {
          mexPrintf(" %ld", line_cnt);
          fflush(NULL);
        }
      }
    }
  }

  fclose(fid);  

  if( mxIsSparse( data_X) )
    jcs[line_cnt] = k;

/*  mexPrintf("\n");*/

  if(verb)
    mexPrintf("Leaving svmlight reading function.\n");

  exitflag = 0;

clean_up:

  mxFree(line);
  mxFree(feat_val);
  mxFree(feat_idx);

  return(exitflag);
}


/*----------------------------------------------------------------------
 Compute area under ROC (1st class label[i]==1; 2nd class label[i] != 1).
  ----------------------------------------------------------------------*/
double compute_auc(double *score, int *label, uint32_t nData)
{
  double *sorted_score = NULL;
  double *sorted_lab = NULL;
  uint32_t i;
  uint32_t neg, pos;
  double auc = -1;

  sorted_score = mxCalloc(nData, sizeof(double));
  if( sorted_score == NULL ) {
      mexPrintf("Not enough memmory to allocate sorted_score when computing AUC.");
      goto clean_up;
  }

  sorted_lab = mxCalloc(nData, sizeof(double));
  if( sorted_lab == NULL )
  {
      mexPrintf("Not enough memmory to allocate sorted_lab when computing AUC.");
      goto clean_up;
  }

  for(i=0; i < nData; i++)
    if(label[i] == 1) sorted_lab[i] = 1.0; else sorted_lab[i] = 0.0;


  memcpy(sorted_score,score,sizeof(double)*nData);

  qsort_data(sorted_score, sorted_lab, nData);

  pos = 0;
  neg = 0;
  auc = 0;

  for(i = 0; i < nData; i++)
  {
    if(sorted_lab[i] ==1.0 )
    {
      pos ++;
    }
    else
    {
      neg ++;
      auc += (double)pos;
    }
  }
  auc = 1 - auc/((double)neg*(double)pos);

clean_up:
  mxFree(sorted_score);
  mxFree(sorted_lab);

  return(auc);
}

