/*=================================================================================
 * LBP_FEATURES_sparse computes pyramid of LBP features for each defined window in input images.
 * 
 * Synopsis:
 *   [Feature,CroppedWin]= lbp_features_sparse(Images,imSize,Wins,winSize,height_of_pyramid,verb) 
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
#include "msvc-compat.h"
#include <math.h>
#include <mex.h>

#include "liblbp.h"

#define INDEX(ROW,COL,NUM_ROWS) ((COL)*(NUM_ROWS)+(ROW))

/*======================================================================
  Main code plus interface to Matlab.
========================================================================*/

void mexFunction( int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[] )
{
  uint32_t i;
  double * tmp;
  uint32_t * lbp;
  mxArray * lbparray;
  int verb;

  if(nrhs < 5 || nrhs > 6)
     mexErrMsgTxt("Improper number of input arguments.\n\n"
                  "LBP_FEATURES_SPARSE computes pyramid of LBP features for each defined window in input images.\n\n"
                  "Synopsis: \n"
                  "  [ Features, sparse_lbp ]= lbppyr_features(Images, imSize, Wins, winSize, height_of_pyramid, verb) \n"
                  "\n"
                  "  Input: \n"
                  "    Images [(im_H*im_W) x nImages (uint8)]\n"
                  "    imSize [2 x 1 (double)] imSize = [im_H im_W]\n"
                  "    Wins [4 x nExamples (uint32)]  [image_idx; x1; y1;mirror] (1-based)\n"
                  "    winSize [2 x 1 (double)] [win_H win_W]\n"
                  "    height_of_pyramid [1 x 1 (double)]\n"
                  "    verb [1x1] \n" 
                  "  Output: \n"
                  "    Features [nDims x nExamples (sparse logical)]\n");

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
  uint32_t nDim = liblbp_pyr_get_dim(win_H,win_W,nPyramids);
  uint32_t nDimLBP = liblbp_pyr_get_dim(win_H,win_W,nPyramids)/256;

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
  
  //lbp = (uint32_t*)mxGetPr(mxCreateNumericMatrix(nDimLBP, nData, mxINT32_CLASS, mxREAL));
  lbparray = mxCreateNumericMatrix(nDimLBP, nData, mxINT32_CLASS, mxREAL);
  lbp = (uint32_t*)mxGetPr(lbparray);

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
              win[cnt0++] = img_ptr[INDEX(y,x,im_H)];
            }

        }
        else
        {
          for(x=x1+win_W-1; x >= x1; x--)
            for(y=y1; y < y1+win_H; y++)
            {
              win[cnt0++] = img_ptr[INDEX(y,x,im_H)];
            }
        }
    
        liblbp_pyr_features_sparse(&lbp[nDimLBP*i], nDimLBP, win, win_H, win_W);
  }
 
  mxFree(win);
  
  mwSize nzmax = nDimLBP*nData;
  plhs[0] = mxCreateSparseLogicalMatrix(nDim, nData, nzmax);
  mwIndex *ir, *jc;
  mxLogical *pr = mxGetLogicals(plhs[0]);
  ir = mxGetIr(plhs[0]);
  jc = mxGetJc(plhs[0]);
  for (mwSize i = 0; i < nzmax; ++i)
  {
      ir[i] = lbp[i];
      pr[i] = 1;
  }
  for (mwSize i = 0; i < nData+1; ++i)
  {
      jc[i] = i*nDimLBP;
  }
  
  // return also lbp_sparse_indices or destroy it
  if (nlhs == 2)
  {
      plhs[1] = lbparray;
  } else {
      mxDestroyArray(lbparray);
  }

  return;
}
