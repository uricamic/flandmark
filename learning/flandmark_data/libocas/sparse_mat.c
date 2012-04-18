/*-----------------------------------------------------------
 * sparse_mat.c: Reimplementation of Matlab functions for 
 *  dealing with sparse matrices and full matrices (support for 
 *  full matrices added later).
 *
  -----------------------------------------------------------*/

#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <stdlib.h>

#include "sparse_mat.h"


INDEX_TYPE_T mxGetM(const mxArray *array_ptr)
{
  return( array_ptr->m );
}

INDEX_TYPE_T mxGetN(const mxArray *array_ptr)
{
  return( array_ptr->n);
}

void mxDestroyArray(mxArray *array_ptr)
{
  if( array_ptr != NULL)
  {
    free(array_ptr->pr);
    free(array_ptr->ir);
    free(array_ptr->jc);
    free(array_ptr);
  }
}

INDEX_TYPE_T mxGetNZMAX(const mxArray *array_ptr)
{
  return(array_ptr->nzmax);
}

mxArray *mxCreateDoubleMatrix(INDEX_TYPE_T m, INDEX_TYPE_T n, mxComplexity ComplexFlag)
{
  mxArray *array = NULL;

  array = (mxArray*)calloc(1,sizeof(mxArray));
  if(array == NULL)
    return(NULL);

  array->pr = NULL;
  array->ir = NULL;
  array->jc = NULL;
  array->m = m;
  array->n = n;
  array->nzmax = m*n;
  array->sparse = 0;

  array->pr = (double*)calloc(m*n,sizeof(double));
  if(array->pr == NULL)
    return(NULL);

  return(array);
}

mxArray *mxCreateSparse(INDEX_TYPE_T m, INDEX_TYPE_T n, NNZ_TYPE_T nzmax, mxComplexity ComplexFlag)
{
  mxArray *array = NULL;

  array = (mxArray*)calloc(1,sizeof(mxArray));
  if(array == NULL)
    return(NULL);

  array->pr = NULL;
  array->ir = NULL;
  array->jc = NULL;
  array->m = m;
  array->n = n;
  array->nzmax = nzmax;
  array->sparse = 1;

  array->pr = (double*)calloc(nzmax,sizeof(double));
  if(array->pr == NULL)
    goto clean_up;

  array->ir = (INDEX_TYPE_T*)calloc(nzmax,sizeof(INDEX_TYPE_T));
  if(array->ir == NULL)
    goto clean_up;

  array->jc = (INDEX_TYPE_T*)calloc(n+1,sizeof(INDEX_TYPE_T));
  if(array->jc == NULL)
    goto clean_up;

  return(array);

clean_up:
  
  free(array->pr);
  free(array->ir);
  free(array->jc);

  return(NULL);
}

double *mxGetPr(const mxArray *array_ptr)
{
  return(array_ptr->pr);
}

INDEX_TYPE_T *mxGetIr(const mxArray *array_ptr)
{
  return(array_ptr->ir);
}

INDEX_TYPE_T *mxGetJc(const mxArray *array_ptr)
{
  return(array_ptr->jc);
}

int mxIsSparse(const mxArray *array_ptr)
{
  return(array_ptr->sparse);
}
