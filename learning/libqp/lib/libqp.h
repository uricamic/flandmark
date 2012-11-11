/*-----------------------------------------------------------------------
 * libqp.h: Library for Quadratic Programming optimization.
 *
 * The library provides two solvers:
 *   1. Solver for QP task with simplex constraints.
 *      See function ./lib/libqp_splx.c for definition of the QP task. 
 *
 *   2. Solver for QP task with box constraints and a single linear 
 *      equality constraint. 
 *      See function ./lib/libqp_gsmo.c for definiton of the QP task. 
 *  
 * Copyright (C) 2006-2008 Vojtech Franc, xfrancv@cmp.felk.cvut.cz
 * Center for Machine Perception, CTU FEL Prague
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public 
 * License as published by the Free Software Foundation; 
 * Version 3, 29 June 2007
 *-------------------------------------------------------------------- */

#ifndef libqp_h
#define libqp_h

#include <stdint.h>
#include <math.h>

#ifdef LIBQP_MATLAB
#include "mex.h"
#define LIBQP_PLUS_INF mxGetInf()
#define LIBQP_CALLOC(x,y) mxCalloc(x,y)
#define LIBQP_FREE(x) mxFree(x)
#else
#define LIBQP_PLUS_INF (-log(0.0))
#define LIBQP_CALLOC(x,y) calloc(x,y)
#define LIBQP_FREE(x) free(x)
#endif

#define LIBQP_INDEX(ROW,COL,NUM_ROWS) ((COL)*(NUM_ROWS)+(ROW))
#define LIBQP_MIN(A,B) ((A) > (B) ? (B) : (A))
#define LIBQP_MAX(A,B) ((A) < (B) ? (B) : (A))
#define LIBQP_ABS(A) ((A) < 0 ? -(A) : (A))

/* QP solver return value */
typedef struct {
  uint32_t nIter;       /* number of iterations */ 
  double QP;            /* primal objective value */ 
  double QD;            /* dual objective value */  
  int8_t exitflag;      /* -1 ... not enough memory 
                            0 ... nIter >= MaxIter 
                            1 ... QP - QD <= TolRel*ABS(QP)
                            2 ... QP - QD <= TolAbs
                            3 ... QP <= QP_TH
                            4 ... eps-KKT conditions satisfied */
} libqp_state_T; 

/* QP solver for tasks with simplex constraints */
libqp_state_T libqp_splx_solver(const double* (*get_col)(uint32_t),
                  double *diag_H,
                  double *f,
                  double *b,
                  uint32_t *I,
                  uint8_t *S,
                  double *x,
                  uint32_t n,
                  uint32_t MaxIter,
                  double TolAbs,
                  double TolRel,
                  double QP_TH,
                  void (*print_state)(libqp_state_T state));

/* Generalized SMO algorithm */
libqp_state_T libqp_gsmo_solver(const double* (*get_col)(uint32_t),
            double *diag_H,
            double *f,
            double *a,
            double b,
            double *LB,
            double *UB,
            double *x,
            uint32_t n,
            uint32_t MaxIter,
            double TolKKT,
            void (*print_state)(libqp_state_T state));


#endif /* libqp_h */
