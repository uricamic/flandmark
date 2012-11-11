/*********************************************************************************
 *
 *  QPSPLX is a standalone program which solves a convex Quadratic Programming 
 *   task with simplex constraints.
 *
 *  Usage: qpsplx input_file output_file
 *
 *  The QP task to solve is defined as
 *
 *    min QP(x) := 0.5*x'*H*x + f'*x  
 *     x
 * 
 * subject to:   
 *   sum_{i in I_k} x[i] == b[k]  for all k such that S[k] == 0 
 *   sum_{i in I_k} x[i] <= b[k]  for all k such that S[k] == 1
 *                             x(i) >= 0 for all i=1:n
 *   
 *  where I_k = { i | I[i] == k}, k={1,...,m}.
 *      
 *  INPUT_FILE is a text file of a simple format which describes:
 *     H [n x n] Symmetric positive semi-definite matrix.
 *     f [n x 1] Vector.
 *     b [m x 1] Vector of positive values.
 *     I [n x 1] Vector of indices containing integers 1,2,...,K (K <= n).
 *     S [m x 1] Vector defining constraint sign: 0.. equality, 1..inequality.
 *  
 *   OUTPUT_FILE is a text file which contains:
 *     QP [1 x 1] Primal objective value.
 *     QD [1 x 1] Dual objective value.
 *     nIter [1 x 1] Number of iterations.
 *     exitflag [1 x 1] Indicates which stopping condition has been used:
 *       -1  ... Not enough memory.
 *       0  ... Maximal number of iteations reached: nIter >= MaxIter.
 *       1  ... Relarive tolerance reached: QP-QD <= abs(QP)*TolRel
 *       2  ... Absolute tolerance reached: QP-QD <= TolAbs
 *       3  ... Objective value reached threshold: QP <= QP_TH.
 *      x [n x 1] Solution vector.
 *
 *
 *   See mtf.h for a description of the used text file format.
 *
 *********************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

#include "mtf.h"
#include "../lib/libqp.h"

double *mat_H = NULL;
uint32_t nVars;

double get_time()
{
	struct timeval tv;
	if (gettimeofday(&tv, NULL)==0)
		return tv.tv_sec+((double)(tv.tv_usec))/1e6;
	else
		return 0.0;
}

const double *get_col_of_mat_H( uint32_t i )
{
  return( &mat_H[ nVars*i ] );
}

int main(int argc, char *argv[])
{
  FILE *input_fid = NULL;
  FILE *output_fid = NULL;
  char *input_fname;  
  char *output_fname; 
  int nMatrices;
  int i, j, idx;
  mtf_matrix_t *matrix_array = NULL;
  uint32_t nConstr;
  double *diag_H = NULL;
  double *vec_f = NULL;
  double *vec_b = NULL;
  uint32_t *vec_I = NULL;
  uint8_t *vec_S = NULL;
  double *vec_x = NULL;
  libqp_state_T exitflag;
  uint32_t MaxIter = 0xFFFFFFFF;
  double TolAbs = 0.0;                  
  double TolRel = 1e-9;
  double QP_TH = -LIBQP_PLUS_INF;
  double load_time = 0;
  double solver_time = 0;
  double total_time = get_time();

  if(argc != 3)
  {
    fprintf(stderr,"Usage: qpsplx input_file output_file\n");
    return(0);
  }
    
  input_fname = argv[1];
  output_fname = argv[2];

  fprintf(stdout,"Input file: %s\n", input_fname);
  fprintf(stdout,"Output file: %s\n", output_fname);

  input_fid = fopen(input_fname, "r");
  if(input_fid == NULL) 
  {
    perror("fopen error: ");
    fprintf(stderr,"Cannot open input file.");
    return(0);
  }

  load_time = get_time();
  nMatrices = mtf_read_file(input_fid, &matrix_array);
  load_time = get_time()-load_time;

  if(nMatrices < 1) 
  {
    fprintf(stderr,"Error when reading matrices from file.\n");
    mtf_print_error( nMatrices );
    goto cleanup;
  }

  if(nMatrices != 5)
  {
    fprintf(stderr,"Input file must contain five matrices/vectors.");
    goto cleanup_free_matrices;
  }

  /*  fprintf(stdout,"Number of loaded matrices: %d\n", nMatrices);
  for(i=0; i < nMatrices; i++ )
  {    
    fprintf(stdout,"name: %s  nRows: %d  nCols: %d\n", 
            matrix_array[i].name, matrix_array[i].nRows, matrix_array[i].nCols);
  }
  */

  /*--------------------------------------------------*/
  /* Get input matrices and vectors                   */
  /*--------------------------------------------------*/
  idx = mtf_get_index("H",matrix_array,nMatrices);
  if(idx < 0 )
  {
    fprintf(stderr,"Matrix H not found.\n");
    goto cleanup_free_matrices;    
  }
  mat_H = matrix_array[idx].value;
  nVars = matrix_array[idx].nCols;

  idx = mtf_get_index("f",matrix_array,nMatrices);
  if(idx < 0 )
  {
    fprintf(stderr,"Vector f not found.\n");
    goto cleanup_free_matrices;    
  }
  vec_f = matrix_array[idx].value;
  
  idx = mtf_get_index("b",matrix_array,nMatrices);
  if(idx < 0 )
  {
    fprintf(stderr,"Vector b not found.\n");
    goto cleanup_free_matrices;    
  }
  vec_b = matrix_array[idx].value;      
  nConstr = LIBQP_MAX(matrix_array[idx].nCols,matrix_array[idx].nRows);
  
  idx = mtf_get_index("I",matrix_array,nMatrices);
  if(idx < 0 )
  {
    fprintf(stderr,"Vector I not found.\n");
    goto cleanup_free_matrices;    
  }
  vec_I = calloc(nVars,sizeof(uint32_t));
  if(vec_I == NULL) {
    fprintf(stderr,"Not enough memory.\n");
    goto cleanup_free_matrices;
  }
  for(i=0; i < nVars; i++)
    vec_I[i] = (uint32_t)matrix_array[idx].value[i];


  idx = mtf_get_index("S",matrix_array,nMatrices);
  if(idx < 0 )
  {
    fprintf(stderr,"Vector S not found.\n");
    goto cleanup_free_matrices;    
  }
  vec_S = calloc(nConstr,sizeof(uint8_t));      
  if(vec_S == NULL) {
    fprintf(stderr,"Not enough memory.\n");
    goto cleanup_free_matrices;
  }
  for(i=0; i < nConstr; i++)
    vec_S[i] = (uint8_t)matrix_array[idx].value[i];


  diag_H = calloc(nVars,sizeof(double));
  if(diag_H == NULL) {
    fprintf(stderr,"Not enough memory.\n");
    goto cleanup_free_matrices;
  }

  for(i=0; i < nVars; i++)
    diag_H[i] = mat_H[LIBQP_INDEX(i,i,nVars)];


  /*--------------------------------------------------*/
  /* Create initial feasible solution                 */
  /*--------------------------------------------------*/
  vec_x = calloc(nVars, sizeof(double)); 
  if( vec_x == NULL )
  {
    fprintf(stderr,"Not enough memory.\n");
    goto cleanup_free_matrices;
  }
  for(i=0; i < nVars; i++)
    vec_x[i] = 0.0;

  for(i=0; i < nConstr; i++ )
  {
    if(vec_S[i] == 0)
    {
      for(j=0; j < nVars; j++)
      {
        if(vec_I[j] == i+1)
        {
          vec_x[j] = vec_b[i];
          break;
        }
      }
    }
  }

  printf("nVar: %d  nConstr: %d\n", nVars, nConstr);

  /*--------------------------------------------------*/
  /* Call the QP solver and print results             */
  /*--------------------------------------------------*/

  printf("Solving QP..."); fflush(stdout); 
  solver_time = get_time();

  exitflag = libqp_splx_solver(get_col_of_mat_H, diag_H, vec_f, vec_b,
                               vec_I, vec_S, vec_x, nVars, MaxIter,
                               TolAbs, TolRel, QP_TH, NULL);

  solver_time = get_time()-solver_time;
  printf("done.\n");

  total_time= get_time() - total_time;

  printf("QP: %f\n", exitflag.QP);
  printf("QD: %f\n", exitflag.QD);
  printf("nIter: %d\n", exitflag.nIter);
  printf("exitflag: %d\n", exitflag.exitflag);
  printf("load_time: %f\n", load_time);
  printf("solver_time: %f\n", solver_time);
  printf("total_time: %f\n", total_time);

  /*--------------------------------------------------*/
  /* Write result to file                             */
  /*--------------------------------------------------*/

  output_fid = fopen(output_fname, "w+");
  if(output_fid == NULL) 
  {
    perror("fopen error: ");
    fprintf(stderr,"Cannot open input file.");
    goto cleanup_free_matrices;
  }

  fprintf(output_fid,"6\n");
  fprintf(output_fid,"QP 1 1\n");
  fprintf(output_fid,"%.12f\n",exitflag.QP);

  fprintf(output_fid,"QD 1 1\n");
  fprintf(output_fid,"%.12f\n",exitflag.QD);

  fprintf(output_fid,"nIter 1 1\n");
  fprintf(output_fid,"%d\n",exitflag.nIter);

  fprintf(output_fid,"exitflag 1 1\n");
  fprintf(output_fid,"%d\n",exitflag.exitflag);

  fprintf(output_fid,"solver_time 1 1\n");
  fprintf(output_fid,"%.12f\n",solver_time);

  fprintf(output_fid,"x %d 1\n", nVars);
  for(i=0; i < nVars; i++ )
    fprintf(output_fid,"%.12f\n", vec_x[i]);

  fclose(output_fid); 

cleanup_free_matrices: 

  free(diag_H);
  free(vec_I); 
  free(vec_S); 
  free(vec_x); 

  for(i=0; i < nMatrices; i++ )  
    free(matrix_array[i].value); 
  
cleanup:  
  fclose(input_fid); 

  return(0); 
}

