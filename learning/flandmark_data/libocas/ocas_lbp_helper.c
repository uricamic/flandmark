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

#include "libocas.h"
#include "ocas_lbp_helper.h"
#include "liblbp.h"

uint8_t *Images;
uint32_t nImages;
uint32_t win_H;
uint32_t win_W;
uint32_t im_H;
uint32_t im_W;
uint32_t nPyramids;
uint32_t *croped_window;
uint32_t *Wins;
uint32_t win_H, win_W;

uint32_t nDim, nData;
double *data_y;
/*double *full_A;*/
int64_t *full_A;
double *W;
double *oldW;
/*double *new_a;*/
/*int32_t *new_a;*/
int64_t *new_a;

double *A0;
double W0;
double oldW0;
double X0;


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


/*-------------------------------------------------------------------------
  sq_norm_W = full_update_W( t ) does the following:

  W = oldW*(1-t) + t*W;
  sq_norm_W = W'*W;
---------------------------------------------------------------------------*/
double full_update_W( double t, void* user_data )
{
  uint32_t j;
  double sq_norm_W;         

/*  mexPrintf("double full_update_W()\n");*/

  W0 = oldW0*(1-t) + t*W0;
  sq_norm_W = W0*W0;

  for(j=0; j <nDim; j++) {
    W[j] = oldW[j]*(1-t) + t*W[j];
    sq_norm_W += W[j]*W[j];
  }          

  return( sq_norm_W );
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
  uint32_t x1, y1, idx, x, y, cnt, mirror;
  uint8_t *img_ptr;

  /*  ptr = mxGetPr(data_X);*/

/*  memset(new_a, 0, sizeof(double)*nDim);*/
  memset(new_a, 0, sizeof(new_a[0])*nDim);

  for(i=0; i < cut_length; i++)
  {

    idx = Wins[LIBOCAS_INDEX(0,new_cut[i],4)]-1;
    x1  = Wins[LIBOCAS_INDEX(1,new_cut[i],4)]-1;
    y1  = Wins[LIBOCAS_INDEX(2,new_cut[i],4)]-1;
    mirror = Wins[LIBOCAS_INDEX(3,new_cut[i],4)];

    img_ptr = &Images[idx*im_H*im_W];
    
    cnt=0;
    if(mirror==0)
    {
      for(x=x1; x < x1+win_W; x++)
        for(y=y1; y < y1+win_H; y++)
          croped_window[cnt++] = img_ptr[LIBOCAS_INDEX(y,x,im_H)];
    }
    else
    {
      for(x=x1+win_W-1; x >= x1; x--)
        for(y=y1; y < y1+win_H; y++)
          croped_window[cnt++] = img_ptr[LIBOCAS_INDEX(y,x,im_H)];
    }
    
    if(data_y[new_cut[i]] == +1) {
/*      lbppyr_addvec(new_a,croped_window);*/
      liblbp_pyr_addvec(new_a,nDim,croped_window,win_H,win_W);
    }
    else
    {
/*      lbppyr_subvec(new_a,croped_window);*/
      liblbp_pyr_subvec(new_a,nDim,croped_window,win_H,win_W);
    }
      
    A0[nSel] += X0*data_y[new_cut[i]];
  }

  /*
  for(i=0; i < cut_length; i++) {
    for(j=0; j < nDim; j++ ) {
      new_a[j] += ptr[LIBOCAS_INDEX(j,new_cut[i],nDim)];
    }

    A0[nSel] += X0*data_y[new_cut[i]];    
  }
  */

  /* compute new_a'*new_a and insert new_a to the last column of full_A */
  sq_norm_a = A0[nSel]*A0[nSel];
  for(j=0; j < nDim; j++ ) 
  {
    sq_norm_a += (double)(new_a[j]*new_a[j]);
    full_A[LIBOCAS_INDEX(j,nSel,nDim)] = new_a[j];
  }

  new_col_H[nSel] = sq_norm_a;
  for(i=0; i < nSel; i++) {
    double tmp = A0[nSel]*A0[i];

    for(j=0; j < nDim; j++ ) {
      tmp += (double)(new_a[j]*full_A[LIBOCAS_INDEX(j,i,nDim)]);
    }
    new_col_H[i] = tmp;
  }

/*  mxFree( new_a );*/

  return 0;
}

/*----------------------------------------------------------------------
  full_compute_output( output ) does the follwing:

  output = data_X'*W;
  ----------------------------------------------------------------------*/
int full_compute_output( double *output, void* user_data )
{
  uint32_t i, j, x1, y1, idx, cnt, x,y, mirror;
  double tmp;
  uint8_t *img_ptr;

/*  mexPrintf("full_compute_output()\n");*/

  /*  ptr = mxGetPr( data_X );*/

  for(i=0; i < nData; i++) 
  { 
    tmp = data_y[i]*X0*W0;

    idx = Wins[LIBOCAS_INDEX(0,i,4)]-1;
    x1  = Wins[LIBOCAS_INDEX(1,i,4)]-1;
    y1  = Wins[LIBOCAS_INDEX(2,i,4)]-1;
    mirror  = Wins[LIBOCAS_INDEX(3,i,4)];

    img_ptr = &Images[idx*im_H*im_W];
    
    cnt=0;
    if(mirror==0)
    {
      for(x=x1; x < x1+win_W; x++)
        for(y=y1; y < y1+win_H; y++)
          croped_window[cnt++] = img_ptr[LIBOCAS_INDEX(y,x,im_H)];
    }
    else
    {
      for(x=x1+win_W-1; x >= x1; x--)
        for(y=y1; y < y1+win_H; y++)
          croped_window[cnt++] = img_ptr[LIBOCAS_INDEX(y,x,im_H)];
    }
    
/*    output[i] = data_y[i]*(X0*W0 + lbppyr_dotprod(W,croped_window));*/
    output[i] = data_y[i]*(X0*W0 + liblbp_pyr_dotprod(W,nDim,croped_window,win_H,win_W));

    /*
    for(j=0; j < nDim; j++ ) {
      tmp += W[j]*ptr[LIBOCAS_INDEX(j,i,nDim)];
    }
    output[i] = tmp;
    */
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
        W[j] += alpha[i]*(double)(full_A[LIBOCAS_INDEX(j,i,nDim)]);
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

