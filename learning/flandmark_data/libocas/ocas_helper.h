/*-----------------------------------------------------------------------
 * ocas_helper.h: Implementation of helper functions for the OCAS solver.
 *
 * It supports both sparse and dense matrices and loading data from
 * the SVM^light format.
 *-------------------------------------------------------------------- */

#ifndef _ocas_helper_h
#define _ocas_helper_h

#include <stdint.h>

#ifdef LIBOCAS_MATLAB

#include <mex.h>

#if !defined(MX_API_VER) || MX_API_VER<0x07040000
#define mwSize int
#define INDEX_TYPE_T int
#define mwIndex int
#else
#define INDEX_TYPE_T mwSize
#endif

#else

#define mwSize int
#define mwIndex int

#include "sparse_mat.h"

#endif

typedef struct {
  double **value;
  uint32_t **index;
  uint32_t *nz_dims;    
} cutting_plane_buf_T;


extern mxArray *data_X;
extern uint32_t nDim, nData, nY;
extern double *data_y;
extern cutting_plane_buf_T sparse_A;
extern double *full_A;
extern double *W;
extern double *oldW;
extern double *new_a;

extern double *A0;
extern double W0;
extern double oldW0;
extern double X0;

double get_time(void);
void ocas_print(ocas_return_value_T value);
void ocas_print_null(ocas_return_value_T value);

void mul_sparse_col(double alpha, mxArray *sparse_mat, uint32_t col);
void add_sparse_col(double *full_vec, mxArray *sparse_mat, uint32_t col);
void subtract_sparse_col(double *full_vec, mxArray *sparse_mat, uint32_t col);
double dp_sparse_col(double *full_vec, mxArray *sparse_mat, uint32_t col);
double sparse_update_W( double t, void* user_data );
int sparse_add_new_cut( double *new_col_H, 
                         uint32_t *new_cut, 
                         uint32_t cut_length, 
                         uint32_t nSel, 
                         void* user_data );
int sparse_compute_output( double *output, void* user_data );
void sparse_compute_W( double *sq_norm_W, double *dp_WoldW, double *alpha, uint32_t nSel, void* user_data );

void msvm_sparse_compute_W( double *sq_norm_W, double *dp_WoldW, double *alpha, uint32_t nSel, void* user_data );
int msvm_sparse_add_new_cut( double *new_col_H, uint32_t *new_cut, uint32_t nSel, void* user_data );
int msvm_sparse_compute_output( double *output, void* user_data );

double msvm_full_update_W( double t, void* user_data );
int msvm_full_add_new_cut( double *new_col_H, uint32_t *new_cut, uint32_t nSel, void* user_data);
int msvm_full_compute_output( double *output, void* user_data );
void msvm_full_compute_W( double *sq_norm_W, double *dp_WoldW, double *alpha, uint32_t nSel, void* user_data );


double full_update_W( double t, void* user_data );
int full_compute_output( double *output, void* user_data );
void full_compute_W( double *sq_norm_W, double *dp_WoldW, double *alpha, uint32_t nSel, void* user_data );
int full_add_new_cut( double *new_col_H, 
                       uint32_t *new_cut, 
                       uint32_t cut_length, 
                       uint32_t nSel,
                       void* user_data);

int load_svmlight_file(char *fname, int verb);
double compute_auc(double *score, int *label, uint32_t nData);
int load_regconsts(char *fname, double **vec_C, uint32_t *len_vec_C, int verb);

int qsort_data(double* value, double* data, uint32_t size);

void destroy_parallel_ocas(void);
int init_parallel_ocas(int number_of_threads);
int parallel_sparse_compute_output( double *output, void* user_data );
int parallel_qsort_data(double* value, double* data, uint32_t size);
int parallel_sparse_add_new_cut( double *new_col_H,uint32_t *new_cut, uint32_t cut_length,uint32_t nSel,void* user_data );


#endif
