/*=================================================================
 * LBPPYR computes LBP features on scale-pyramid of input image.
 * 
 * Synopsis:
 *  F = lbppyr_mat(I, P)
 * where
 *  I [H x W (uint8)] is input image.
 *  P [1 x 1 (double)] is height of the scale-pyramid.
 *  F [N x 1 (uint8)] LBP features stucked to a column vector. 
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

  uint8_t *I, K;
  uint32_t *P;
  uint32_t *J;
  int w,h,x,y,i,j,N;
  int ww,hh; 
  uint32_t center;
  uint32_t offset;
  uint8_t pattern;

  if( nrhs != 2 )
    mexErrMsgTxt("Two input arguments required.\n\n"
                 "LBPPYR computes LBP features on scale-pyramid of input image.\n"
                 "Synopsis:\n"
                 "  F = lbppyr(I, P)\n"
                 "where\n"
                 "  I [H x W (uint8)] is input image.\n"
                 "  P [1 x 1 (double)] is height of the scale-pyramid.\n"
                 "  F [N x 1 (uint8)] LBP features stucked to a column vector. \n");

  I = (uint8_t*)mxGetPr(prhs[0]);
  K = (uint8_t)mxGetScalar(prhs[1]);

  h = mxGetM(prhs[0]);
  w = mxGetN(prhs[0]);

  /*  printf("h: %d\n", h);
  printf("w: %d\n", w);
  printf("K: %d\n", K);
  */

  /* count number of LBPs */
  for(ww=w, hh=h, N=0, i=0; i < K && MIN(ww,hh) >= 3; i++)
  {
    N += (ww-2)*(hh-2);

    if(ww % 2) ww--;
    if(hh % 2) hh--;
    ww = ww/2;
    hh = hh/2;
  }
  K = i;
  N = 256*N;
  
/*  printf("N: %d (=256*%d)\n", N,N/256);*/
/*  printf("K: %d\n",K);*/
  /*  printf("%d/%d = %d*%d + %d \n", b,a, a, b/a, b % a);*/
  
  plhs[0] = mxCreateNumericMatrix(N, 1, mxUINT32_CLASS, mxREAL);
  P = (uint32_t*)mxGetPr(plhs[0]);

  J = mxCalloc(h*w, sizeof(uint32_t));
  if(J ==NULL)
    mexErrMsgTxt("Not enough memory.");

  for(x=0; x < w; x++)
    for(y=0; y < h; y++)
      J[INDEX(y,x,h)] = I[INDEX(y,x,h)];
  
  for(ww=w, hh=h, i=0, offset = 0; i < K; i++)
  {
    for(x=1; x < ww-1; x++)
    {
      for(y=1; y< hh-1; y++)
      {
        pattern = 0;
        center = J[INDEX(y,x,h)];
        if(J[INDEX(y-1,x-1,h)] < center) pattern = pattern | 0x01;
        if(J[INDEX(y-1,x,h)] < center)   pattern = pattern | 0x02;
        if(J[INDEX(y-1,x+1,h)] < center) pattern = pattern | 0x04;
        if(J[INDEX(y,x-1,h)] < center)   pattern = pattern | 0x08;
        if(J[INDEX(y,x+1,h)] < center)   pattern = pattern | 0x10;
        if(J[INDEX(y+1,x-1,h)] < center) pattern = pattern | 0x20;
        if(J[INDEX(y+1,x,h)] < center)   pattern = pattern | 0x40;
        if(J[INDEX(y+1,x+1,h)] < center) pattern = pattern | 0x80;

        /*        P[cnt++] = pattern; */
        P[offset+pattern] = 1;
        offset += 256;
      }
    }

    if( i < K-1 )
    {
      if(ww % 2 == 1) ww--;
      if(hh % 2 == 1) hh--;

      ww = ww/2;
      for(x=0; x < ww; x++)
        for(j=0; j < hh; j++)
          J[INDEX(j,x,h)] = J[INDEX(j,2*x,h)] + J[INDEX(j,2*x+1,h)];

      hh = hh/2;
      for(y=0; y < hh; y++)
        for(j=0; j < ww; j++)
          J[INDEX(y,j,h)] = J[INDEX(2*y,j,h)] + J[INDEX(2*y+1,j,h)];
    }
    
  }
  
  return;
}

