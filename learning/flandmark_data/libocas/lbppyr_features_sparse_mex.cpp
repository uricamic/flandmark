/*=================================================================================
 * LBPPYR_FEATURES computes pyramid of LBP features for each defined window in input images.
 * 
 * Synopsis:
 *   [Feature,CroppedWin]= lbppyr_features(Images,imSize,Wins,winSize,height_of_pyramid,verb) 
 * Input: 
 *   Images [(im_H*im_W) x nImages (uint8)]
 *   imSize [2 x 1 (double)] imSize = [im_H im_W]
 *   Wins [4 x nExamples (uint32)]  [image_index; top_left_col; top_left_row; mirror]
 *   winSize [2 x 1 (double)] [win_H win_W]
 *   height_of_pyramid [1 x 1 (double)] 
 * Output: 
 *   Features [nDims x nExamples (sparse logical)]
 *   CroppedWin [(win_H*win_W) x nExamples]
 *
 *======================================================================================*/ 

#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <mex.h>
#include "liblbp.h"

#define INDEX(ROW,COL,NUM_ROWS) ((COL)*(NUM_ROWS)+(ROW))


/*======================================================================
  Main code plus interface to Matlab.
========================================================================*/

void mexFunction( int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[] )
{
  uint32_t i;//, j;
  double *tmp;
  char *TmpFeatures;
  uint8_t *CroppedWin;
  int verb;

  if(nrhs < 5 || nrhs > 6)
     mexErrMsgTxt("Improper number of input arguments.\n\n"
                  "LBPPYR_FEATURES computes pyramid of LBP features for each defined window in input images.\n\n"
                  "Synopsis: \n"
                  "  [Features,CroppedWin]= lbppyr_features(Images,imSize,Wins,winSize,height_of_pyramid,verb) \n"
                  "\n"
                  "  Input: \n"
                  "    Images [(im_H*im_W) x nImages (uint8)]\n"
                  "    imSize [2 x 1 (double)] imSize = [im_H im_W]\n"
                  "    Wins [4 x nExamples (uint32)]  [image_idx; x1; y1;mirror] (1-based)\n"
                  "    winSize [2 x 1 (double)] [win_H win_W]\n"
                  "    height_of_pyramid [1 x 1 (double)]\n"
                  "    verb [1x1] \n" 
                  "  Output: \n"
                  "    Features [nDims x nExamples (sparse logical)]\n"
                  "    CroppedWin [(win_H*win_W) x nExamples]\n");

  if(nrhs == 6)
    verb = (int)mxGetScalar(prhs[5]);
  else
    verb = 1;

  uint8_t * Images = (uint8_t*)mxGetPr(prhs[0]);
  uint32_t nImages = mxGetN(prhs[0]);

  tmp = (double*)mxGetPr(prhs[1]);
  uint32_t im_H = (uint32_t)tmp[0];
  uint32_t im_W = (uint32_t)tmp[1];

  if(mxGetM(prhs[0]) != im_H*im_W)
    mexErrMsgTxt("Dimension of Images does not match to im_H*im_W.");

  uint32_t *Wins = (uint32_t*)mxGetPr(prhs[2]);

  tmp = (double*)mxGetPr(prhs[3]);
  uint16_t win_H = (uint16_t)tmp[0];
  uint16_t win_W = (uint16_t)tmp[1];

  uint16_t nPyramids = (uint32_t)mxGetScalar(prhs[4]);
  uint32_t nDim =  liblbp_pyr_get_dim(win_H,win_W,nPyramids);

  uint32_t nData = mxGetN(prhs[2]);

  if(verb)
  {
    mexPrintf("Input data:\n"
              "   # of images     : %d\n"
              "   image height    : %d\n"
              "   image width     : %d\n",
              nImages, im_H, im_W);

    mexPrintf("Feature represenation:\n"
              "   base window height        : %d\n"
              "   base window width         : %d\n"
              "   nPyramids                 : %d\n"
            "   # of virtual examples     : %d\n"
              "   # of features per example : %d\n",
              win_H, win_W, nPyramids, nData, nDim);
  }

  /* learned weight vector */
/*  plhs[0] = mxCreateNumericMatrix(nDim,nData, mxINT64_CLASS, mxREAL);*/
/*  Features = (int64_t*)mxGetPr(plhs[0]);*/
  
  TmpFeatures = (char*)mxGetPr(mxCreateNumericMatrix(nDim, nData, mxINT8_CLASS, mxREAL));

  plhs[1] = mxCreateNumericMatrix(win_H*win_W, nData, mxUINT8_CLASS, mxREAL);
  CroppedWin = (uint8_t*)mxGetPr(plhs[1]);


  uint32_t cnt, cnt0, mirror,x,x1,y,y1,idx;
  uint32_t *win;
  uint8_t *img_ptr;
  
  win = (uint32_t*)mxCalloc(win_H*win_W,sizeof(uint32_t));
  if(win == NULL) 
    mexErrMsgTxt("Not enough memory for croped_window.");

  cnt=0;
  for(i=0; i < nData; i++)
  {

    idx = Wins[INDEX(0,i,4)]-1;
    x1  = Wins[INDEX(1,i,4)]-1;
    y1  = Wins[INDEX(2,i,4)]-1;
    mirror = Wins[INDEX(3,i,4)];

    img_ptr = &Images[idx*im_H*im_W];
 
    cnt0 = 0;
    if(mirror==0)
    {
      for(x=x1; x < x1+win_W; x++)
        for(y=y1; y < y1+win_H; y++)
        {
          CroppedWin[cnt++] = img_ptr[INDEX(y,x,im_H)];
          win[cnt0++] = img_ptr[INDEX(y,x,im_H)];
        }
         
    }
    else
    {
      for(x=x1+win_W-1; x >= x1; x--)
        for(y=y1; y < y1+win_H; y++)
        {
          CroppedWin[cnt++] = img_ptr[INDEX(y,x,im_H)];
          win[cnt0++] = img_ptr[INDEX(y,x,im_H)];
        }
    }
    
/*    lbppyr_features(uint32_t *vec, uint32_t *win);*/
/*    liblbp_pyr_features( &Features[nDim*i], win);*/
/*    mexPrintf("i=%d, nDim=%d, im_H=%d, im_W=%d\n", i,nDim,im_H,im_W);*/
    liblbp_pyr_features(&TmpFeatures[nDim*i], nDim, win, win_H, win_W);
  }
 
  int cntOnesColumn = 0;
  for (int i = 0; i < nDim; ++i)
  {
	  cntOnesColumn += TmpFeatures[i];
  }
 
  mwIndex *irs, *jcs;
  mxLogical *sr;
  mwSize nzmax = (mwSize)(nData*cntOnesColumn);
  plhs[0] = mxCreateSparseLogicalMatrix(nDim, nData, nzmax);
  sr = mxGetLogicals(plhs[0]);
  irs = mxGetIr(plhs[0]);
  jcs = mxGetJc(plhs[0]);

  // fill in sparse logical matrix
  mwIndex k = 0;
  mwIndex jidx;
  
  for (jidx = 0; jidx < nData; ++jidx)
  {
      mwSize i;
      jcs[jidx] = k;
      for (i = 0; i < nDim; ++i)
      {
          // copy non-zero element
          if (TmpFeatures[i])
          {
              //mexPrintf("Non zero element found. Adding to sparse matrix...\n");
              sr[k] = 1;
              irs[k] = i;
              ++k;
          }
      }
      TmpFeatures += nDim;
  }
  jcs[nData] = k;
  
  
  mxFree(win);

  return;
}

