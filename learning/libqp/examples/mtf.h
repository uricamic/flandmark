/********************************************************************************

 A library for loading dense matrices stored in a text file format.

 Example text file:           description (not contained in the file)
 +--------------------------+ 
 | 2                        | Number of matrices stored in the file.
 | A 2 3                    | NAME NUM_ROWS NUM_COLS
 | 1                        | Elements of the matrix stored column-wise
 | 2                        | 
 | 3                        |
 | 4                        |
 | 5                        |
 | 6                        |
 | b 2 1                    | NAME NUM_ROWS NUM_COLS
 | 10                       | Elements of the matrix stored column-wise
 | 20                       |
 ---------------------------+

 NAME is a string
 NUM_ROWS is a positive integer
 NUM_COLS is a positive integer
 Elements are doubles.

 The example file (insided the frame) describes two matrices 

   A = [1 3 5]
       [2 4 6]

   b = [10]
       [20]

 Copyright (C) 2006-2008 Vojtech Franc, xfrancv@cmp.felk.cvut.cz
 Center for Machine Perception, CTU FEL Prague

********************************************************************************/

#ifndef mtf_h
#define mtf_h

#include <stdint.h>

#define MTF_PREMATURE_EOF -1
#define MTF_NOT_ENOUGH_MEMORY -2
#define MTF_DEFECTIVE_FORMAT -3

#define MTF_MAX_LINE_LENGTH 1024

typedef struct 
{ 
  char name[MTF_MAX_LINE_LENGTH];
  uint32_t nRows;
  uint32_t nCols;
  double *value;
  
} mtf_matrix_t;

void mtf_print_error( int error_code );

int mtf_read_file(FILE *fid, mtf_matrix_t **matrix_array );

int mtf_get_index(char *item, mtf_matrix_t *matrix_array, int nMatrices);


#endif 
