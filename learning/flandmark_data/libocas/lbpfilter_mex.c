/*=================================================================
 * LBPFILTER computes LBP for each 3x3 subwindow in input image.
 *
 * Synopsis:
 *  F = lbpfilter(I)
 * where
 *  I [h x w (double)] is input image.
 *  F [h x w (double)] is image of LBP responses;
 *   the border colums and rows, for which LBP is not defined, are set 0.
 *
 *=================================================================*/

#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <mex.h>
#include <time.h>
#include <errno.h>

#define MIN(A,B) ((A) > (B) ? (B) : (A))
#define MAX(A,B) ((A) < (B) ? (B) : (A))
#define ABS(A) ((A) < 0 ? -(A) : (A))
#define INDEX(ROW,COL,NUM_ROWS) ((COL)*(NUM_ROWS)+(ROW))


/*======================================================================
  Main code plus interface to Matlab.
========================================================================*/

void mexFunction( int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[] )
{

  double *I, *P;
  int w,h, x,y;
  double center;
  unsigned char pattern;

  if( nrhs != 1 )
    mexErrMsgTxt("One input argument required.\n\n"
                 "LBPFILTER computes LBP for each 3x3 subwindow in input image.\n"
                 "\n"
                 "Synopsis:\n"
                 "  F = lbpfilter(I)\n"
                 "where\n"
                 "  I [h x w (double)] is input image.\n"
                 "  F [h x w (double)] is image of LBP responses; \n"
                 "    the border colums and rows, for which LBP is not defined, are set 0.\n"
                 );

  I = (double*)mxGetPr(prhs[0]);
  h = mxGetM(prhs[0]);
  w = mxGetN(prhs[0]);

  plhs[0] = mxCreateDoubleMatrix(h,w,mxREAL);
  P = (double*)mxGetPr(plhs[0]);

  for(x=1; x < w-1; x++)
  {
    for(y=1; y< h-1; y++)
    {
      pattern = 0;
      center = I[INDEX(y,x,h)];
      if(I[INDEX(y-1,x-1,h)] < center) pattern = pattern | 0x01;
      if(I[INDEX(y-1,x,h)] < center)   pattern = pattern | 0x02;
      if(I[INDEX(y-1,x+1,h)] < center) pattern = pattern | 0x04;
      if(I[INDEX(y,x-1,h)] < center)   pattern = pattern | 0x08;
      if(I[INDEX(y,x+1,h)] < center)   pattern = pattern | 0x10;
      if(I[INDEX(y+1,x-1,h)] < center) pattern = pattern | 0x20;
      if(I[INDEX(y+1,x,h)] < center)   pattern = pattern | 0x40;
      if(I[INDEX(y+1,x+1,h)] < center) pattern = pattern | 0x80;

      P[INDEX(y,x,h)] = (double)pattern; 
    }
  }
  
  return;
}

