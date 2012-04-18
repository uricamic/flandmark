/*=================================================================
 * COMPUTE_ERRORS computes classification error and area under ROC.
 *
 *  Synopsis:
 *    [error,auc] = compute_errors(score,true_labels)
 *
 *  Input:
 *   scores [1 x n] scores of binary classifier; the classifier decides 
 *      according to the sign of the score.
 *   true_labels [1 x n] true labels; 1st class label = 1, 
 *                                    2nd class label != 1
 *  Output:
 *   error [1x1] number of mis-classifications
 *   auc [1x1] Area under ROC 
 *
 *=================================================================*/

#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <mex.h>
#include <sys/time.h>
#include <time.h>
#include <errno.h>

#include "libocas.h"
#include "ocas_helper.h"

#if !defined(MX_API_VER) || MX_API_VER<0x07040000
#define mwSize int
#define mwIndex int
#endif

/*======================================================================
  Main code.
========================================================================*/

void mexFunction( int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[] )
{
  int *true_labels_int;
  uint32_t n1, n2, i;
  double auc, error;

  double *scores, *true_labels;

  if( nrhs != 2 )
    mexErrMsgTxt("Two input arguments must be passed.\n"
                 "\n"
                 "[error,auc] = compute_errors(score,true_labels)\n"
                 "Inputs: score [N x 1 (double)] positive class score >= 0 otherwise negative class \n"
                 "        true_labels [N x 1 (double)] positive class 1; negative class ~= 1\n"
                 "Outputs: error [1x1] number of misclassifications\n"
                 "         auc [1x1] Arrea Under ROC\n");

  scores = mxGetPr(prhs[0]);
  true_labels = mxGetPr(prhs[1]);

  n1 = LIBOCAS_MAX(mxGetM(prhs[0]),mxGetN(prhs[0]));
  n2 = LIBOCAS_MAX(mxGetM(prhs[1]),mxGetN(prhs[1]));

  if(n1 != n2)
    mexErrMsgTxt("The input vectors must be of the same size.");


  true_labels_int = mxCalloc(n1,sizeof(int));
  if(true_labels_int == NULL)
    mexErrMsgTxt("Not enough memory.");

  error = 0;
  for(i = 0; i < n1; i++)
  {
    true_labels_int[i] = (int)true_labels[i];

    if((true_labels[i] == 1 && scores[i] <= 0) || (true_labels[i] != 1 && scores[i] > 0))
      error++;
  }

  auc = compute_auc(scores, true_labels_int, n1);

  mxFree(true_labels_int);

  plhs[0] = mxCreateDoubleScalar(error);
  plhs[1] = mxCreateDoubleScalar(auc);

  return;
}

