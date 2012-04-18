#include "flandmark_detector.h"

#include <mex.h>
#include <ctime>

// MATLAB entry point
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs < 2 || nrhs > 3)
    {
        mexErrMsgTxt("Improper number of input arguments.\n\n"
                "FLANDMARK_DETECTOR detects facial landmarks.\n"
                "Synopsis: \n"
                "flandmark_detector(face_image, fname)\n"
                "\n"
                "Input: \n"
                "  face_image [im_H * im_W (uint8)] column vector \n"
                "  fname      [string]\n"
                "\n"
             );
    }
    
    if (!mxIsUint8(prhs[0]))
    {
        mexErrMsgTxt("Improper argument type \"face_image\"\n\n"
                "  face image must be [im_H * im_W] uint8 matrix.\n"
              );
    }
    
    uint8_t * face_image = (uint8_t*)mxGetPr(prhs[0]);
    char * fname = mxArrayToString(prhs[1]);
    //clock_t start, stop;
    double start, ms;

    /* check input */
    const mwSize *face_dim;
    face_dim = mxGetDimensions(prhs[0]);
    if (face_image == NULL)
    {
        mexErrMsgTxt("Improper argument type \"face_image\"\n\n"
                "  face image must be [im_H x im_W] uint8 matrix.\n"
              );
    }
    if (fname == NULL)
    {
        mexErrMsgTxt("Improper argument type \"fname\"\n\n"
                "  fname must be a string.\n"
              );
    }
    
    /* read model */
    start = (double)cvGetTickCount();
    FLANDMARK_Model * model = flandmark_init(fname);
    if (model == 0)
    {
        mexErrMsgTxt("Error opening file with model structure.\n");
    }
    start = (double)cvGetTickCount() - start;
    ms = cvRound( start / ((double)cvGetTickFrequency() * 1000.0) );
    
    if (face_dim[0] != model->data.imSize[0]*model->data.imSize[1])
    {
        mexErrMsgTxt("Improper argument type \"face_image\"\n\n"
                "  Dimensions of face_image does not match dimensions in model dat file.\n"
              );
    }
    
    mexPrintf("Structure model loaded in %.6f ms\n", ms);

    /* TEST - call detector */
    start = (double)cvGetTickCount();
    float * landmarks = (float*)malloc(2*model->data.options.M*sizeof(float));
    flandmark_detect_base(face_image, model, landmarks);
    start = (double)cvGetTickCount() - start;
    ms = cvRound( start / ((double)cvGetTickFrequency() * 1000.0) );
    
    mexPrintf("Landmarks positions estimated in %.6f ms\n", ms);

    start = (double)cvGetTickCount();
    plhs[0] = mxCreateNumericMatrix(2, model->data.options.M, mxDOUBLE_CLASS, mxREAL);
    double * temp = (double*)mxGetPr(plhs[0]);
    for (int i = 0; i < 2*model->data.options.M; ++i)
    {
        temp[i] = (double)landmarks[i];
    }
    start = (double)cvGetTickCount() - start;
    ms = cvRound( start / ((double)cvGetTickFrequency() * 1000.0) );
    
    free(landmarks);

    mexPrintf("Output for MATLAB computed in %.6f ms\n", ms);

    /* cleanup */
    mexPrintf("Cleaning structure model... ");
    flandmark_free(model);
    mexPrintf(" done.\n");
}
