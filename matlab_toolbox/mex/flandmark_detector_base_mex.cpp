/*
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * Written (W) 2012 Michal Uricar
 * Copyright (C) 2012 Michal Uricar
 */

#include "flandmark_detector.h"

#include <mex.h>
#include <ctime>

// MATLAB entry point
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs < 2 || nrhs > 3)
    {
        mexErrMsgTxt("Improper number of input arguments.\n\n"
                "FLANDMARK_DETECTOR_BASE detects facial landmarks.\n"
                "Synopsis: \n"
                "flandmark_detector(face_image, fname)\n"
                "\n"
                "Input: \n"
                "  face_image [im_H * im_W (uint8)] column vector \n"
                "  model [(uint8)] data from flandmark_model.dat \n"
                "Output: \n"
                "  landmarks [2 x M] matrix (DOUBLE) x, y coordinates of landmarks in original image\n"
                "\n"
             );
    }
    
    if (!mxIsUint8(prhs[0]))
    {
        mexErrMsgTxt("Improper argument type \"face_image\"\n\n"
                "  face image must be [im_H * im_W] uint8 matrix.\n"
              );
    }
    
    if (!mxIsUint8(prhs[1]))
    {
        mexErrMsgTxt("Improper argument type \"face_image\"\n\n"
                "  face image must be uint8 array loaded by flandmark_load_model.\n"
              );
    }

    uint8_t * face_image = (uint8_t*)mxGetPr(prhs[0]);

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

    FLANDMARK_Model * model = (FLANDMARK_Model*)mxGetPr(prhs[1]);

    const mwSize *mx_dim = mxGetDimensions(prhs[0]);
    if (mx_dim[0]*mx_dim[1] != model->data.options.bw[0]*model->data.options.bw[1])
    {
        mexErrMsgTxt("Improper argument type \"bbox\"\n\n"
                     "   bbox must be [1 x 4 (int)] int32 vector.\n"
            );
    }

    start = (double)cvGetTickCount();

    plhs[0] = mxCreateNumericMatrix(2, model->data.options.M, mxDOUBLE_CLASS, mxREAL);
    double *landmarks = (double*)mxGetPr(plhs[0]);
    flandmark_detect_base(face_image, model, landmarks);

    start = (double)cvGetTickCount() - start;
    ms = cvRound( start / ((double)cvGetTickFrequency() * 1000.0) );

    mexPrintf("Landmarks positions estimated in %.6f ms\n", ms);
}
