/*-----------------------------------------------------------------------
 * libocas.h: Implementation of the OCAS solver for training 
 *            linear SVM classifiers.
 *  
 * Copyright (C) 2008, 2009 Vojtech Franc, xfrancv@cmp.felk.cvut.cz
 *                          Soeren Sonnenburg, soeren.sonnenburg@first.fraunhofer.de
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public 
 * License as published by the Free Software Foundation; 
 *-------------------------------------------------------------------- */


#ifndef libocas_h
#define libocas_h

#include <stdint.h>
#include "libqp.h"

#ifdef LIBOCAS_MATLAB

#include "mex.h"
#define LIBQP_MATLAB
#define LIBOCAS_PLUS_INF mxGetInf()
#define LIBOCAS_CALLOC(x,y) mxCalloc(x,y)
#define LIBOCAS_FREE(x) mxFree(x)

#else

#define LIBOCAS_PLUS_INF (-log(0.0))
#define LIBOCAS_CALLOC(x,y) calloc(x,y)
#define LIBOCAS_FREE(x) free(x)

#endif

#define LIBOCAS_INDEX(ROW,COL,NUM_ROWS) ((COL)*(NUM_ROWS)+(ROW))
#define LIBOCAS_MIN(A,B) ((A) > (B) ? (B) : (A))
#define LIBOCAS_MAX(A,B) ((A) < (B) ? (B) : (A))
#define LIBOCAS_ABS(A) ((A) < 0 ? -(A) : (A))



typedef struct {
  uint32_t nIter;        /* number of iterations */
  uint32_t nCutPlanes;   /* number of cutitng buffered planes */
  uint32_t nNZAlpha;     /* number of non-zero Lagrangeans (effective number of CPs) */
  uint32_t trn_err;      /* number of training errors */
  double Q_P;            /* primal objective value */
  double Q_D;            /* dual objective value */
  double output_time;    /* time spent in computing outputs */
  double sort_time;      /* time spent in sorting */
  double add_time;       /* time spent in adding examples to compute cutting planes */
  double w_time;         /* time spent in computing parameter vector  */
  double qp_solver_time; /* time spent in inner QP solver  */
  double ocas_time;      /* total time spent in svm_ocas_solver */
  double print_time;     /* time spent in ocas_print function */
  int8_t qp_exitflag;    /* exitflag from the last call of the inner QP solver */
  int8_t exitflag;       /*  1 .. ocas.Q_P - ocas.Q_D <= TolRel*ABS(ocas.Q_P) 
                             2 .. ocas.Q_P - ocas.Q_D <= TolAbs 
                             3 .. ocas.Q_P <= QPBound
                             4 .. optimization time >= MaxTime 
                            -1 .. ocas.nCutPlanes >= BufSize 
                            -2 .. not enough memory for the solver */
} ocas_return_value_T;

/* binary linear SVM solver */
ocas_return_value_T svm_ocas_solver(
         double C,            /* regularizarion constant */
         uint32_t nData,      /* number of exmaples */
         double TolRel,       /* halts if 1-Q_P/Q_D <= TolRel */
         double TolAbs,       /* halts if Q_P-Q_D <= TolRel */
         double QPBound,      /* halts if QP <= QPBound */
         double MaxTime,      /* maximal time in seconds spent in optmization */
         uint32_t BufSize,    /* maximal number of buffered cutting planes  */
         uint8_t Method,      /* 0..standard CP (SVM-Perf,BMRM), 1..OCAS */
         void (*compute_W)(double*, double*, double*, uint32_t, void*),
         double (*update_W)(double, void*),
         int (*add_new_cut)(double*, uint32_t*, uint32_t, uint32_t, void*),
         int (*compute_output)( double*, void* ),
         int (*sort)(double*, double*, uint32_t),
         void (*ocas_print)(ocas_return_value_T),
         void* user_data);

/* binary linear SVM solver which allows using different C for each example*/
ocas_return_value_T svm_ocas_solver_difC(
         double *C,           /* regularizarion constants for each example */
         uint32_t nData,      /* number of exmaples */
         double TolRel,       /* halts if 1-Q_P/Q_D <= TolRel */
         double TolAbs,       /* halts if Q_P-Q_D <= TolRel */
         double QPBound,      /* halts if QP <= QPBound */
         double MaxTime,      /* maximal time in seconds spent in optmization */
         uint32_t BufSize,    /* maximal number of buffered cutting planes  */
         uint8_t Method,      /* 0..standard CP (SVM-Perf,BMRM), 1..OCAS */
         void (*compute_W)(double*, double*, double*, uint32_t, void*),
         double (*update_W)(double, void*),
         int (*add_new_cut)(double*, uint32_t*, uint32_t, uint32_t, void*),
         int (*compute_output)( double*, void* ),
         int (*sort)(double*, double*, uint32_t),
         void (*ocas_print)(ocas_return_value_T),
         void* user_data);

/* multi-class (Singer-Crammer formulation) linear SVM solver */
ocas_return_value_T msvm_ocas_solver(
            double C,
            double *data_y,
            uint32_t nY,
            uint32_t nData, 
            double TolRel,
            double TolAbs,
            double QPBound,
            double MaxTime,
            uint32_t _BufSize,
            uint8_t Method,
            void (*compute_W)(double*, double*, double*, uint32_t, void*),
            double (*update_W)(double, void*),
            int (*add_new_cut)(double*, uint32_t*, uint32_t, void*),
            int (*compute_output)(double*, void* ),
            int (*sort)(double*, double*, uint32_t),
			void (*ocas_print)(ocas_return_value_T),
			void* user_data);

#endif /* libocas_h */

