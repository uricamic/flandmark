/*
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * Written (W) 2012 Michal Uricar
 * Copyright (C) 2012 Michal Uricar
 */

#include <highgui.h>
#include <cstring>
#include <cmath>

#include "flandmark_detector.h"

int main( int argc, char** argv ) 
{
    char flandmark_window[] = "flandmark_simple_example";
    double t;
    int ms;
    int * bbox = (int*)malloc(4*sizeof(int));
    
    if (argc < 6)
    {
      fprintf(stderr, "Usage: flandmark_1 <path_to_input_image> <face_bbox - 4int> [<path_to_output_image>]\n");
      exit(1);
    }
    
    //cvNamedWindow(flandmark_window, 0 );
    
    t = (double)cvGetTickCount();
    FLANDMARK_Model * model = flandmark_init("flandmark_model.dat");
    if (model == 0)
    {
        printf("Structure model wasn't created. Corrupted file flandmark_model.dat?\n");
        exit(1);
    }
    t = (double)cvGetTickCount() - t;
    ms = cvRound( t / ((double)cvGetTickFrequency() * 1000.0) );
    printf("Structure model loaded in %d ms.\n", ms);
    
    
    // input image
    IplImage *img = cvLoadImage(argv[1]);
    if (img == NULL)
    {
      //fprintf(stderr, "Wrong path to image. Exiting...\n");
      fprintf(stderr, "Cannot open image %s. Exiting...\n", argv[1]);
      exit(1);
    }

    // convert image to grayscale
    IplImage *img_grayscale = cvCreateImage(cvSize(img->width, img->height), IPL_DEPTH_8U, 1);
    cvCvtColor(img, img_grayscale, CV_BGR2GRAY);
    
    // face bbox
    bbox[0] = ::atoi(argv[2]);
    bbox[1] = ::atoi(argv[3]);
    bbox[2] = ::atoi(argv[4]);
    bbox[3] = ::atoi(argv[5]);


    // call flandmark_detect
    t = (double)cvGetTickCount();
    double * landmarks = (double*)malloc(2*model->data.options.M*sizeof(double));
    if(flandmark_detect(img_grayscale, bbox, model, landmarks))
    {
        printf("Error during detection.\n");
    }
    t = (double)cvGetTickCount() - t;
    ms = cvRound( t / ((double)cvGetTickFrequency() * 1000.0) );
    printf("Landmarks detected in %d ms.\n", ms);

    cvRectangle(img, cvPoint(bbox[0], bbox[1]), cvPoint(bbox[2], bbox[3]), CV_RGB(255,0,0) );
    cvRectangle(img, cvPoint(model->bb[0], model->bb[1]), cvPoint(model->bb[2], model->bb[3]), CV_RGB(0,0,255) );
    cvCircle(img, cvPoint((int)landmarks[0], (int)landmarks[1]), 3, CV_RGB(0, 0,255), CV_FILLED);
    for (int i = 2; i < 2*model->data.options.M; i += 2)
    {
        cvCircle(img, cvPoint(int(landmarks[i]), int(landmarks[i+1])), 3, CV_RGB(255,0,0), CV_FILLED);

    }
    printf("detection = \t[");
    for (int ii = 0; ii < 2*model->data.options.M; ii+=2)
    {
            printf("%.2f ", landmarks[ii]);
    }
    printf("]\n");
    printf("\t\t[");
    for (int ii = 1; ii < 2*model->data.options.M; ii+=2)
    {
            printf("%.2f ", landmarks[ii]);
    }
    printf("]\n");

    cvShowImage(flandmark_window, img);
    cvWaitKey(0);


    if (argc == 3)
    {
      printf("Saving image to file %s...\n", argv[2]);
      cvSaveImage(argv[2], img);
    }
    
    // cleanup
    cvDestroyWindow(flandmark_window);
    cvReleaseImage(&img);
    cvReleaseImage(&img_grayscale);
    free(landmarks);
    free(bbox);
    flandmark_free(model);
}
