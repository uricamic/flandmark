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
				"getNormalizedFrame creates normalized image frame given the image and detected face bounding box.\n"
				"Synopsis: \n"
				"[face_image, bb] = getNormalizedFrame(image, bbox, model.data.options)\n"
				"\n"
				"Input: \n"
				"  image [h * w (uint8)] column vector \n"
				"  bbox  [1 x 4 (double)]\n"
				"  model.data.options  [1 x 1 (struct)]\n"
				"Output: \n"
				"  face_image [bw[0] X bw[1] (uint8)] matrix\n"
				"  bb [1 x 4 (double)] vector\n"
				"\n"
				);
	}

	const mwSize *dims;
	double *bw = 0, *bw_margin = 0, *p_bbox = 0;
	int bbox[4];
	IplImage *ipl_image;
	FLANDMARK_Model *model = (FLANDMARK_Model*)malloc(sizeof(FLANDMARK_Model));

	// check input
	if (!mxIsStruct(prhs[2]))
	{
		mexErrMsgTxt("Improper argument type \"options\".\n\n"
				"  model must be a structure containing bw and bw_margin.\n"
				);
	}
	mxArray *in_bw = mxGetField(prhs[2], 0, "bw");
	mxArray *in_bw_margin = mxGetField(prhs[2], 0, "bw_margin");
	if (in_bw == NULL || in_bw_margin == NULL)
	{
		mexErrMsgTxt("Improper argument type \"options\".\n\n"
				"  model must be a structure containing bw and bw_margin.\n"
				);
	}

	bw = (double*)mxGetPr(in_bw);
	dims = mxGetDimensions(in_bw);
	if (dims[0] != 2 || dims[1] != 1)
	{
		mexErrMsgTxt("Improper argument type \"options\".\n\n"
				"  options must contain field bw - [2 x 1 (double)] vector.\n"
				);
	}
	model->data.options.bw[0] = bw[0];
	model->data.options.bw[1] = bw[1];

	bw_margin = (double*)mxGetPr(in_bw_margin);
	dims = mxGetDimensions(in_bw);
	if (dims[0] != 2 || dims[1] != 1)
	{
		mexErrMsgTxt("Improper argument type \"options\".\n\n"
				"  options must contain field bw_margin - [2 x 1 (double)] vector.\n"
				);
	}
	model->data.options.bw_margin[0] = bw_margin[0];
	model->data.options.bw_margin[1] = bw_margin[1];

	p_bbox = (double*)mxGetPr(prhs[1]);
	dims = mxGetDimensions(prhs[1]);
	if (p_bbox == 0 || dims[0] != 1 || dims[1] != 4)
	{
		mexErrMsgTxt("Improper argument type \"bbox\".\n\n"
				"  bbox must be vector [1 x 4 (double)] vector.\n"
				);
	}
	bbox[0] = p_bbox[0]; bbox[1] = p_bbox[1];
	bbox[2] = p_bbox[2]; bbox[3] = p_bbox[3];

	if ((int)mxGetNumberOfDimensions(prhs[0]) != 2)
	{
		mexErrMsgTxt("Improper argument type \"image\".\n\n"
				"  image must be uint8 2D matrix (grayscale image).\n"
				);
	}
	dims = mxGetDimensions(prhs[0]);

	ipl_image = mx_to_IplImage(prhs[0]);

	plhs[0] = mxCreateNumericMatrix(bw[1], bw[0], mxUINT8_CLASS, mxREAL);
	uint8_t *p_face_image = (uint8_t*)mxGetPr(plhs[0]);

	plhs[1] = mxCreateNumericMatrix(1, 4, mxDOUBLE_CLASS, mxREAL);
	double *p_bb = (double*)mxGetPr(plhs[1]);

    model->sf = (float*)malloc(2*sizeof(float));

    if (flandmark_get_normalized_image_frame(ipl_image, bbox, p_bb, p_face_image, model))
	{
		//    mexErrMsgTxt("Error ?\n");
        mexErrMsgTxt("Error obtaining normalized image frame!");
	}

	// destroy model (or to be more precise just what was created)
	cvReleaseImage(&ipl_image);
    free(model->sf);
}
