/********************************************************************************

 A library for loading dense matrices stored in a text file format.

 Copyright (C) 2006-2008 Vojtech Franc, xfrancv@cmp.felk.cvut.cz
 Center for Machine Perception, CTU FEL Prague

********************************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mtf.h"

int mtf_get_index(char *item, mtf_matrix_t *matrix_array, int nMatrices)
{
  int i;

  for(i=0; i < nMatrices; i++ )
  {
    if(strcmp(item,matrix_array[i].name) == 0)
      return(i);
  }
  
  return(-1);
}


void mtf_print_error( int error_code )
{
  switch(error_code)
  {
     case MTF_PREMATURE_EOF:
       fprintf(stderr,"Premature end of the file.\n");
       break;

     case MTF_NOT_ENOUGH_MEMORY:
       fprintf(stderr,"Not enought memory.\n");
       break;

     case MTF_DEFECTIVE_FORMAT:
       fprintf(stderr,"Defective format.\n");       
       break;
  }
}

int mtf_read_file(FILE *fid, mtf_matrix_t **matrix_array )
{
  char line[MTF_MAX_LINE_LENGTH];
  char name[MTF_MAX_LINE_LENGTH];
  uint32_t nCols, nRows;
  int num_items_read;
  int exitflag;
  uint32_t i, j;
  int nMatrices;
  double value;

  /* now continue scanning until you reach the end-of-comments */
  do 
  {
    if (fgets(line,MTF_MAX_LINE_LENGTH,fid) == NULL) 
      return MTF_PREMATURE_EOF;
  } 
  while(line[0] == '#' || line[0] == '\n');

  if( sscanf(line, "%d\n", &nMatrices) != 1)
    return MTF_DEFECTIVE_FORMAT;
  
  if(nMatrices <= 0 )
    return MTF_DEFECTIVE_FORMAT;

  *matrix_array = calloc(nMatrices, sizeof(mtf_matrix_t));
  if( matrix_array == NULL)
    return MTF_NOT_ENOUGH_MEMORY;

  for(i=0; i < nMatrices; i++ )
    (*matrix_array)[i].value = NULL;
  
  for(i=0; i < nMatrices; i++ )
  {
      num_items_read = fscanf(fid, "%s %d %d\n", name, &nRows, &nCols); 
      if (num_items_read == EOF) 
      {
        exitflag = MTF_PREMATURE_EOF;
        goto cleanup;
      }
      if( num_items_read != 3 || nRows <=0 || nCols <= 0)
      {
        exitflag = MTF_DEFECTIVE_FORMAT;
        goto cleanup;
      }

      strcpy((*matrix_array)[i].name,name);
      (*matrix_array)[i].nCols = nCols;
      (*matrix_array)[i].nRows = nRows;
      (*matrix_array)[i].value = calloc(nRows*nCols, sizeof(double));
      if( (*matrix_array)[i].value == NULL ) 
      {
        exitflag = MTF_NOT_ENOUGH_MEMORY;
        goto cleanup;
      }

      for(j=0; j < nRows*nCols; j++ )
      {        
        num_items_read = fscanf(fid, "%lg\n", &value); 
        if (num_items_read == EOF) 
        {
          exitflag = MTF_DEFECTIVE_FORMAT;
          goto cleanup;
        }
        if( num_items_read != 1)
        {
          exitflag = MTF_DEFECTIVE_FORMAT;
          goto cleanup;
        }

        (*matrix_array)[i].value[j] = value;
      }
    
  }

  return(nMatrices);

cleanup:

  for(i=0; i < nMatrices; i++ )
    free((*matrix_array)[i].value);
  
  return(exitflag);
}
