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

extern uint8_t *Images;
extern uint32_t nImages;
extern uint32_t win_H;
extern uint32_t win_W;
extern uint32_t im_H;
extern uint32_t im_W;
extern uint32_t nPyramids;
extern uint32_t *croped_window;
extern uint32_t *Wins;

extern uint32_t nDim, nData;
extern double *data_y;
/*extern double *full_A;*/
/*extern int32_t *full_A;*/
extern int64_t *full_A;
extern double *W;
extern double *oldW;
/*extern double *new_a;*/
/*extern int32_t *new_a;*/
extern int64_t *new_a;

extern double *A0;
extern double W0;
extern double oldW0;
extern double X0;

uint32_t lbppyr_get_dim(uint32_t win_H, uint32_t win_W, uint32_t nPyramids);

double get_time(void);
void ocas_print(ocas_return_value_T value);
void ocas_print_null(ocas_return_value_T value);

double full_update_W( double t, void* user_data );
int full_compute_output( double *output, void* user_data );
void full_compute_W( double *sq_norm_W, double *dp_WoldW, double *alpha, uint32_t nSel, void* user_data );
int full_add_new_cut( double *new_col_H, 
                       uint32_t *new_cut, 
                       uint32_t cut_length, 
                       uint32_t nSel,
                       void* user_data);

double compute_auc(double *score, int *label, uint32_t nData);

int qsort_data(double* value, double* data, uint32_t size);

void lbppyr_features(char *vec, uint32_t *win);

#endif
