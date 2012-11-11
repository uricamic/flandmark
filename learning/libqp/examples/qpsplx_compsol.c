/***************************************************************************************

  This program compares two solutions produced by qpsplx solver.

  Usage: qpsplx_compsol solution_file1 solution_file2
 
 ***************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

#include "mtf.h"
#include "../lib/libqp.h"

int main(int argc, char *argv[])
{
  FILE *input_fid1 = NULL;
  FILE *input_fid2 = NULL;
  char *input_fname1;  
  char *input_fname2; 
  int nMatrices1, nMatrices2;
  int i, j, idx;
  mtf_matrix_t *matrix_array1 = NULL;
  mtf_matrix_t *matrix_array2 = NULL;
  double *vec_x1, *vec_x2;
  double QP1, QP2;
  double QD1, QD2;
  double time1, time2;
  double sol_dif;
  uint32_t nIter1, nIter2;
  uint32_t nVars1, nVars2;
  uint32_t nnz1 = 0, nnz2 = 0;

  if(argc != 3)
  {
    fprintf(stderr,"Usage: qpsplx_compsol solution_file1 solution_file2\n");
    return(0);
  }
    
  input_fname1 = argv[1];
  input_fname2 = argv[2];

  fprintf(stdout,"Solution file 1: %s\n", input_fname1);
  fprintf(stdout,"Solution file 2: %s\n", input_fname2);

  input_fid1 = fopen(input_fname1, "r");
  if(input_fid1 == NULL) 
  {
    perror("fopen error: ");
    fprintf(stderr,"Cannot open solution file 1.");
    goto cleanup;
  }

  input_fid2 = fopen(input_fname2, "r");
  if(input_fid2 == NULL) 
  {
    perror("fopen error: ");
    fprintf(stderr,"Cannot open solution file 2.");
    goto cleanup;
  }


  /*--------------------------------------------------*/
  /* Get solution from file 1                         */
  /*--------------------------------------------------*/

  nMatrices1 = mtf_read_file(input_fid1, &matrix_array1);

  if(nMatrices1 < 1) 
  {
    fprintf(stderr,"Error when reading matrices from solution file 1.\n");
    mtf_print_error( nMatrices1 );
    goto cleanup;
  }

  if(nMatrices1 != 6)
  {
    fprintf(stderr,"Input file 1 must contain six vectors/scalars.");
    goto cleanup_free_matrices;
  }

  idx = mtf_get_index("QP",matrix_array1,nMatrices1);
  if(idx < 0 )
  {
    fprintf(stderr,"Scalar QP not found in file 1.\n");
    goto cleanup_free_matrices;    
  }
  QP1 = matrix_array1[idx].value[0];

  idx = mtf_get_index("QD",matrix_array1,nMatrices1);
  if(idx < 0 )
  {
    fprintf(stderr,"Scalar QD not found in file 1.\n");
    goto cleanup_free_matrices;    
  }
  QD1 = matrix_array1[idx].value[0];

  idx = mtf_get_index("solver_time",matrix_array1,nMatrices1);
  if(idx < 0 )
  {
    fprintf(stderr,"Scalar solver_time not found in file 1.\n");
    goto cleanup_free_matrices;    
  }
  time1 = matrix_array1[idx].value[0];

  idx = mtf_get_index("nIter",matrix_array1,nMatrices1);
  if(idx < 0 )
  {
    fprintf(stderr,"Scalar nIter not found in file 1.\n");
    goto cleanup_free_matrices;    
  }
  nIter1 = (uint32_t)matrix_array1[idx].value[0];

  idx = mtf_get_index("x",matrix_array1,nMatrices1);
  if(idx < 0 )
  {
    fprintf(stderr,"Vector x not found in file 1.\n");
    goto cleanup_free_matrices;    
  }
  vec_x1 = matrix_array1[idx].value;
  nVars1 = LIBQP_MAX(matrix_array1[idx].nCols,matrix_array1[idx].nRows);

  /*--------------------------------------------------*/
  /* Get solution from file 2                         */
  /*--------------------------------------------------*/

  nMatrices2 = mtf_read_file(input_fid2, &matrix_array2);

  if(nMatrices2 < 1) 
  {
    fprintf(stderr,"Error when reading matrices from solution file 2.\n");
    mtf_print_error( nMatrices2 );
    goto cleanup;
  }

  if(nMatrices2 != 6)
  {
    fprintf(stderr,"Input file 2 must contain six vectors/scalars.");
    goto cleanup_free_matrices;
  }

  idx = mtf_get_index("QP",matrix_array2,nMatrices2);
  if(idx < 0 )
  {
    fprintf(stderr,"Scalar QP not found in file 2.\n");
    goto cleanup_free_matrices;    
  }
  QP2 = matrix_array2[idx].value[0];

  idx = mtf_get_index("QD",matrix_array2,nMatrices2);
  if(idx < 0 )
  {
    fprintf(stderr,"Scalar QD not found in file 2.\n");
    goto cleanup_free_matrices;    
  }
  QD2 = matrix_array2[idx].value[0];

  idx = mtf_get_index("solver_time",matrix_array2,nMatrices2);
  if(idx < 0 )
  {
    fprintf(stderr,"Scalar solver_time not found in file 2.\n");
    goto cleanup_free_matrices;    
  }
  time2 = matrix_array2[idx].value[0];

  idx = mtf_get_index("nIter",matrix_array2,nMatrices2);
  if(idx < 0 )
  {
    fprintf(stderr,"Scalar nIter not found in file 2.\n");
    goto cleanup_free_matrices;    
  }
  nIter2 = (uint32_t)matrix_array2[idx].value[0];

  idx = mtf_get_index("x",matrix_array2,nMatrices2);
  if(idx < 0 )
  {
    fprintf(stderr,"Vector x not found in file 2.\n");
    goto cleanup_free_matrices;    
  }
  vec_x2 = matrix_array2[idx].value;
  nVars2 = LIBQP_MAX(matrix_array2[idx].nCols,matrix_array2[idx].nRows);

  /*--------------------------------------------------------------*/
  /* Compare solutions                                            */
  /*--------------------------------------------------------------*/
  if(nVars1 != nVars2)
  {
    fprintf(stderr,"Length of solution vectors must not differ.\n");
    goto cleanup_free_matrices;        
  }

  sol_dif = 0;
  for(i = 0; i < nVars1; i++)
  {
    sol_dif = LIBQP_MAX(sol_dif, LIBQP_ABS(vec_x1[i]-vec_x2[i]));

    if(vec_x1[i]) nnz1++;
    if(vec_x2[i]) nnz2++;
  }

  printf("                               prim_val      time      nnz    iter       dual_gap\n");
  printf("solution file 1         %.10f     %.3f  %7d   %5d   %.10f\n",
            QP1, time1, nnz1, nIter1, QP1-QD1);
  printf("solution file 2         %.10f     %.3f  %7d   %5d   %.10f\n",
            QP2, time2, nnz2, nIter2, QP2-QD2);  
  printf("Difference in solution vectors: max(abs(x1-x2)) = %.10f\n", sol_dif);


cleanup_free_matrices: 

  if(matrix_array1 != NULL )
  {
    for(i=0; i < nMatrices1; i++ )  
      free(matrix_array1[i].value); 
  }

  if(matrix_array2 != NULL )
  {
    for(i=0; i < nMatrices2; i++ )  
      free(matrix_array2[i].value); 
  }
  
cleanup:  
  fclose(input_fid1); 
  fclose(input_fid2); 

  return(0); 
}

