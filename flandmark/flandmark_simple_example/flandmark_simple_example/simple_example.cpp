///*
//* This program is free software; you can redistribute it and/or modify
//* it under the terms of the GNU General Public License as published by
//* the Free Software Foundation; either version 3 of the License, or
//* (at your option) any later version.
//*
//* Written (W) 2012 Michal Uricar
//* Copyright (C) 2012 Michal Uricar
//*/
//
//#include <highgui.h>
//#include <cstring>
//#include <cmath>
//#include <iostream>
//
//#include <fstream>
//
//#include "flandmark_detector.h"
//
//using namespace std;
//
//#include <sstream>
//#include <vector>
//
//std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
//	std::stringstream ss(s);
//	std::string item;
//	while (std::getline(ss, item, delim)) {
//		elems.push_back(item);
//	}
//	return elems;
//}
//
//
//std::vector<std::string> split(const std::string &s, char delim) {
//	std::vector<std::string> elems;
//	split(s, delim, elems);
//	return elems;
//}
//
//
//
//
//int main()
//{
//	char flandmark_window[] = "flandmark_simple_example";
//	double t;
//	int ms;
//	int * bbox = (int*)malloc(4 * sizeof(int));
//
//	//cvNamedWindow(flandmark_window, 0 );
//
//	t = (double)cvGetTickCount();
//	FLANDMARK_Model * model = flandmark_init("flandmark_model.dat");
//	if (model == 0)
//	{
//		printf("Structure model wasn't created. Corrupted file flandmark_model.dat?\n");
//		exit(1);
//	}
//
//	t = (double)cvGetTickCount() - t;
//	ms = cvRound(t / ((double)cvGetTickFrequency() * 1000.0));
//	printf("Structure model loaded in %d ms.\n", ms);
//
//
//	//// input image
//	//char* input_image = new char[500]; 
//	//input_image = "E:/Users/daniyuu/Lab/flandmark/data/Images/Aaron_Eckhart_0001.jpg";
//
//	IplImage *img = cvLoadImage("E:/Users/daniyuu/Lab/flandmark/data/Images/Aaron_Eckhart_0001.jpg");
//	//if (img == NULL)
//	//{
//	//	//fprintf(stderr, "Wrong path to image. Exiting...\n");
//	//	fprintf(stderr, "Cannot open image %s. Exiting...\n", argv[1]);
//	//	exit(1);
//	//}
//
//	// open file ''example.dat'' for reading and writing
//	
//
//	ifstream in("E:/Users/daniyuu/Lab/flandmark/data/Images/Aaron_Eckhart_0001.det");
//
//	string line;
//	getline(in, line); 
//
//	vector<string> x = split(line, ',');
//
//	// convert image to grayscale
//	IplImage *img_grayscale = cvCreateImage(cvSize(img->width, img->height), IPL_DEPTH_8U, 1);
//	cvCvtColor(img, img_grayscale, CV_BGR2GRAY);
//
//	// face bbox
//	for (int i = 0; i < x.size(); ++i){
//		string s = x[i];
//		stringstream ss;
//		ss << s;
//		int temp;
//		ss >> temp;
//		bbox[i] = temp;
//
//		cout << bbox[i] << endl;
//	}
//	
//
//
//	// call flandmark_detect
//	t = (double)cvGetTickCount();
//	double * landmarks = (double*)malloc(2 * model->data.options.M*sizeof(double));
//	if (flandmark_detect(img_grayscale, bbox, model, landmarks))
//	{
//		printf("Error during detection.\n");
//	}
//	t = (double)cvGetTickCount() - t;
//	ms = cvRound(t / ((double)cvGetTickFrequency() * 1000.0));
//	printf("Landmarks detected in %d ms.\n", ms);
//
//	cvRectangle(img, cvPoint(bbox[0], bbox[1]), cvPoint(bbox[2], bbox[3]), CV_RGB(255, 0, 0));
//	cvRectangle(img, cvPoint(model->bb[0], model->bb[1]), cvPoint(model->bb[2], model->bb[3]), CV_RGB(0, 0, 255));
//	cvCircle(img, cvPoint((int)landmarks[0], (int)landmarks[2]), 3, CV_RGB(0, 0, 255), CV_FILLED);
//	for (int i = 2; i < 2 * model->data.options.M; i += 2)
//	{
//		cout << "options == i == " << i << endl;
//		cvCircle(img, cvPoint(int(landmarks[i]), int(landmarks[i + 1])), 3, CV_RGB(255, 0, 0), CV_FILLED);
//
//	}
//	printf("detection = \t[");
//	for (int ii = 0; ii < 2 * model->data.options.M; ii += 2)
//	{
//		printf("%.2f ", landmarks[ii]);
//	}
//	printf("]\n");
//	printf("\t\t[");
//	for (int ii = 1; ii < 2 * model->data.options.M; ii += 2)
//	{
//		printf("%.2f ", landmarks[ii]);
//	}
//	printf("]\n");
//
//	cvShowImage(flandmark_window, img);
//	cvWaitKey(0);
//
//
//	/*if (argc == 3)
//	{
//		printf("Saving image to file %s...\n", argv[2]);
//		cvSaveImage(argv[2], img);
//	}
//
//*/
//
//	
//	cvSaveImage("E:/Users/daniyuu/Lab/flandmark/record/simple_example/out_Aaron_Eckhart_0001.jpg", img);
//	cout << "here " << endl;
//	system("pause");
//
//	// cleanup
//	cvDestroyWindow(flandmark_window);
//	cvReleaseImage(&img);
//	cvReleaseImage(&img_grayscale);
//	free(landmarks);
//	free(bbox);
//	flandmark_free(model);
//}
//
//
