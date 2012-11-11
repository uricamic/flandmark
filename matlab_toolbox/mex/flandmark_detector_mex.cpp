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

IplImage* mx_to_IplImage(const mxArray *mx_img)
{
    const mwSize *dims;
    int nDims;
    uint8_t *pImg;
    IplImage *img = 0;

    // Get dimensions
    nDims = (int)mxGetNumberOfDimensions(mx_img);
    dims = mxGetDimensions(mx_img);
    CvSize img_size = cvSize(dims[1], dims[0]);

    pImg = (uint8_t*)mxGetPr(mx_img);

    img = cvCreateImage(img_size, IPL_DEPTH_8U, 1);
    for (long x = 0; x < dims[1]; ++x)
    {
        for (long y = 0; y < dims[0]; ++y)
        {
            ((uint8_t*)(img->imageData + img->widthStep*y))[x] = pImg[x*dims[0]+y];
        }
    }

    return img;
}

// MATLAB entry point
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs != 3)
    {
        mexErrMsgTxt("Improper number of input arguments.\n\n"
                "FLANDMARK_DETECTOR detects facial landmarks.\n"
                "Synopsis: \n"
                "flandmark_detector(image, bbox, fname)\n"
                "\n"
                "Input: \n"
                "  image [h * w (uint8)] matrix (UINT8, grayscale) input image \n"
                "  bbox  [1 x 4 (int)] row vector (INT32) detected face bounding box \n"
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

    if (!mxIsUint8(prhs[2]))
    {
        mexErrMsgTxt("Improper argument type \"face_image\"\n\n"
                "  face image must be uint8 array loaded by flandmark_load_model.\n"
              );
    }

    FLANDMARK_Model * model = (FLANDMARK_Model*)mxGetPr(prhs[2]);

    double start, ms;

    /* check input */
    const mwSize *mx_dim;
    //face_dim = mxGetDimensions(prhs[0]);

    IplImage *image = mx_to_IplImage(prhs[0]);
    if (image == NULL)
    {
        mexErrMsgTxt("Improper argument type \"face_image\"\n\n"
                "  face image must be [im_H x im_W] uint8 matrix.\n"
              );
    }
    
    mx_dim = mxGetDimensions(prhs[1]);
    if (mx_dim[0]*mx_dim[1] != 4)
    {
        mexErrMsgTxt("Improper argument type \"bbox\"\n\n"
                     "   bbox must be [1 x 4 (int)] int32 vector.\n"
            );
    }

    /* get bbox form input */
    int *bbox = (int *)mxGetPr(prhs[1]);

    /* call detector */

    start = (double)cvGetTickCount();

    plhs[0] = mxCreateNumericMatrix(2, model->data.options.M, mxDOUBLE_CLASS, mxREAL);
    double *landmarks = (double*)mxGetPr(plhs[0]);
    flandmark_detect(image, bbox, model, landmarks);

    start = (double)cvGetTickCount() - start;
    ms = cvRound( start / ((double)cvGetTickFrequency() * 1000.0) );

    mexPrintf("Landmarks positions estimated in %.6f ms\n", ms);

    // cleanup
    cvReleaseImage(&image);
}
