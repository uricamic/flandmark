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
                "  fname [string] path to the flandmark model data file \n"
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
    
    //uint8_t * image = (uint8_t*)mxGetPr(prhs[0]);
    IplImage *image = mx_to_IplImage(prhs[0]);
    char * fname = mxArrayToString(prhs[2]);
    //clock_t start, stop;
    double start, ms;

    /* check input */
    const mwSize *face_dim;
    face_dim = mxGetDimensions(prhs[0]);
    if (image == NULL)
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
    
    /* get bbox form input */
    int *bbox = (int *)mxGetPr(prhs[1]);

    /* read model */
    start = (double)cvGetTickCount();
    FLANDMARK_Model * model = flandmark_init(fname);
    if (model == 0)
    {
        mexErrMsgTxt("Error opening file with model structure.\n");
    }
    //stop = clock();
    start = (double)cvGetTickCount() - start;
    ms = cvRound( start / ((double)cvGetTickFrequency() * 1000.0) );
    mexPrintf("Structure model loaded in %.6f ms\n", ms);

    /* call detector */
    start = (double)cvGetTickCount();
    float * landmarks = (float*)malloc(2*model->data.options.M*sizeof(float));
    flandmark_detect(image, bbox, model, landmarks);
    start = (double)cvGetTickCount() - start;
    ms = cvRound( start / ((double)cvGetTickFrequency() * 1000.0) );
    
    mexPrintf("Landmarks positions estimated in %.6f ms\n", ms);

    /* Create output */
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
