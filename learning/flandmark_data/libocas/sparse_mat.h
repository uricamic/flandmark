/*-----------------------------------------------------------
 * sparse_mat.h: Reimplementation of Matlab functions for 
 *  dealing with sparse matrices.
 *
  -----------------------------------------------------------*/

#ifndef _sparse_mat_h
#define _sparse_mat_h

#define mxCalloc(x...) calloc(x)
#define mxFree(x...) free(x)
#define mexPrintf(x...) printf(x)

#define INDEX_TYPE_T uint32_t
#define NNZ_TYPE_T uint64_t


typedef struct {
  INDEX_TYPE_T *ir;       
  INDEX_TYPE_T *jc;  
  INDEX_TYPE_T m;
  INDEX_TYPE_T n;
  double *pr;
  NNZ_TYPE_T nzmax;
  int sparse; 
  
} mxArray;


typedef enum {
    mxREAL,
    mxCOMPLEX
} mxComplexity;


INDEX_TYPE_T mxGetM(const mxArray *array_ptr);
INDEX_TYPE_T mxGetN(const mxArray *array_ptr);
void mxDestroyArray(mxArray *array_ptr);
mxArray *mxCreateSparse(INDEX_TYPE_T m, INDEX_TYPE_T n, NNZ_TYPE_T nzmax, mxComplexity ComplexFlag);
mxArray *mxCreateDoubleMatrix(INDEX_TYPE_T m, INDEX_TYPE_T n, mxComplexity ComplexFlag);
double *mxGetPr(const mxArray *array_ptr);
INDEX_TYPE_T *mxGetIr(const mxArray *array_ptr);
INDEX_TYPE_T *mxGetJc(const mxArray *array_ptr);
INDEX_TYPE_T mxGetNZMAX(const mxArray *array_ptr);
int mxIsSparse(const mxArray *array_ptr);

#endif
