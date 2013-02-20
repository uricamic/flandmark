/*
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * Written (W) 2012 Michal Uricar
 * Copyright (C) 2012 Michal Uricar
 */

#include <mex.h>
#include "msvc-compat.h"

#include "flandmark_detector.h"

// MATLAB entry point
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs != 1)
    {
        mexErrMsgTxt("Improper number of input arguments.\n\n"
                "FLANDMARK_LOAD_MODEL loads model for flandmark detector to memory.\n\n"
                "Synopsis: \n"
                "model = flandmark_load_model(fname)\n"
                "\n"
                "Input: \n"
                "  fname [string]\n"
                "Output: \n"
                "  model [binary data]\n"
            );
    }

    int address_offset, dynamic_begin, lbp_pointer, disp_pointer;
    int model_size = 0, dynamic_offset = 0;
    double start, ms;
    char * fname = mxArrayToString(prhs[0]);
    int tsize = 0;
    FLANDMARK_PSIG * PsiGi = NULL;

    // check input
    if (fname == NULL)
    {
        mexErrMsgTxt("Improper argument type \"fname\"\n\n"
                     "  fname must be a string.");
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
    mexPrintf("Structure model loaded in %.6f ms\n", ms);

    uint8_t M = model->data.options.M;

    // Compute size of model_size -------------------------------------------
    model_size += sizeof(FLANDMARK_Model);
    model_size += model->W_ROWS*model->W_COLS*sizeof(double);
    model_size += model->data.options.bw[0]*model->data.options.bw[1]*sizeof(uint8_t);
    model_size += 4*sizeof(double);
    model_size += 2*sizeof(float);
    model_size += sizeof(FLANDMARK_LBP*);
    model_size += M*sizeof(FLANDMARK_LBP);
    for (int i = 0; i < M; ++i)
    {
        model_size += model->data.lbp[i].WINS_ROWS*model->data.lbp[i].WINS_COLS*sizeof(uint32_t);
    }
    model_size += 4*M*sizeof(double);
    model_size += sizeof(FLANDMARK_Options);
    model_size += 4*M*sizeof(int);
    for (int i = 0; i < 3; ++i)
    {
        switch (i)
        {
            case 0: PsiGi = model->data.options.PsiGS0;
            break;
            case 1: PsiGi = model->data.options.PsiGS1;
            break;
            case 2: PsiGi = model->data.options.PsiGS2;
            break;
        }

        tsize = model->data.options.PSIG_ROWS[i]*model->data.options.PSIG_COLS[i];

        for (int idx = 0; idx < tsize; ++idx)
        {
            model_size += sizeof(FLANDMARK_PSIG);
            model_size += PsiGi[idx].ROWS*PsiGi[idx].COLS*sizeof(int);
        }
    }

    mexPrintf("model_size = %d\n", model_size);
    //mexErrMsgTxt("Exitting for sure...\n");

    /* Create uint8 buffer for storing the structure into matlab memory */
    plhs[0] = mxCreateNumericMatrix(model_size, 1, mxUINT8_CLASS, mxREAL);
    uint8_t *buffer = (uint8_t*)mxGetPr(plhs[0]);

    // copy data from structure to memory
    // Statically allocated parts -----------------------------------------------

    // FLANDMARK_Model - static part
    address_offset = 0;
    memcpy(&buffer[address_offset], model, sizeof(FLANDMARK_Model));

    // FLANDMARK_Data - static part
    address_offset = (char*)&(model->data) - (char*)model;
    memcpy(&buffer[address_offset], &model->data, sizeof(FLANDMARK_Data));

    // FLANDMARK_Options - static part
    address_offset = (char*)&(model->data.options) - (char*)model;
    memcpy(&buffer[address_offset], &model->data.options, sizeof(FLANDMARK_Options));

    // Dynamically allocated parts ----------------------------------------------
    dynamic_offset = sizeof(FLANDMARK_Model);
    dynamic_begin = dynamic_offset;

    // model->W
    address_offset = dynamic_offset;
    memcpy(&buffer[address_offset], model->W, model->W_COLS*model->W_ROWS*sizeof(double));
    dynamic_offset += model->W_COLS*model->W_ROWS*sizeof(double);
    // update pointer model->W
    *(double**)&buffer[(char*)&(model->W)-(char*)model]=(double*)&buffer[address_offset];

    // model->normalizedImageFrame
    address_offset = dynamic_offset;
    memcpy(&buffer[address_offset], model->normalizedImageFrame, model->data.options.bw[0]*model->data.options.bw[1]*sizeof(uint8_t));
    dynamic_offset += model->data.options.bw[0]*model->data.options.bw[1]*sizeof(uint8_t);
    // update pointer model->normalizedImageFrame
    *(uint8_t**)&buffer[(char*)&(model->normalizedImageFrame)-(char*)model]=(uint8_t*)&buffer[address_offset];

    // model->bb
    address_offset = dynamic_offset;
    memcpy(&buffer[address_offset], model->bb, 4*sizeof(double));
    dynamic_offset += 4*sizeof(double);
    // update pointer model->bb
    *(double**)&buffer[(char*)&(model->bb)-(char*)model]=(double*)&buffer[address_offset];

    // model->sf
    address_offset = dynamic_offset;
    memcpy(&buffer[address_offset], model->sf, 2*sizeof(float));
    dynamic_offset += 2*sizeof(float);
    // update pointer model->sf
    *(float**)&buffer[(char*)&(model->sf)-(char*)model]=(float*)&buffer[address_offset];

    // model->data.mapTable
    address_offset = dynamic_offset;
    memcpy(&buffer[address_offset], model->data.mapTable, M*4*sizeof(int));
    dynamic_offset += M*4*sizeof(int);
    // update pointer model->data.mapTable
    *(int**)&buffer[(char*)&(model->data.mapTable)-(char*)model]=(int*)&buffer[address_offset];

    // model->data.options.S
    address_offset = dynamic_offset;
    memcpy(&buffer[address_offset], model->data.options.S, M*4*sizeof(int));
    dynamic_offset += M*4*sizeof(int);
    // update pointer model->data.options.S
    *(int**)&buffer[(char*)&(model->data.options.S)-(char*)model]=(int*)&buffer[address_offset];

    // model->data.lbp

    // lbp static part
    lbp_pointer = dynamic_offset;
    for (uint8_t i=0; i < M; ++i)
    {
        address_offset = dynamic_offset;
        memcpy(&buffer[address_offset], &model->data.lbp[i], sizeof(FLANDMARK_LBP));
        dynamic_offset += sizeof(FLANDMARK_LBP);
    }
    // Update pointer
    *(FLANDMARK_LBP**)&buffer[(char*)&(model->data.lbp)-(char*)model]=(FLANDMARK_LBP*)&buffer[lbp_pointer];

    // lbp dynamic part
    for (uint8_t i=0; i < M; ++i)
    {
        address_offset = dynamic_offset;
        memcpy(&buffer[address_offset], model->data.lbp[i].wins, model->data.lbp[i].WINS_COLS*model->data.lbp[i].WINS_ROWS*sizeof(uint32_t));
        dynamic_offset += model->data.lbp[i].WINS_COLS*model->data.lbp[i].WINS_ROWS*sizeof(uint32_t);
        // update pointer
        *(uint32_t**)&buffer[lbp_pointer+i*sizeof(FLANDMARK_LBP)+(long long)((char*)&(model->data.lbp[i].wins))-(long long)((char*)&(model->data.lbp[i]))]=(uint32_t*)&buffer[address_offset];
    }

    // model->data.options.PsiGSi

    for (int i=0; i < 3; ++i)
    {
        address_offset = dynamic_offset;
        switch(i)
        {
            case 0: PsiGi = model->data.options.PsiGS0;
                *(FLANDMARK_PSIG**)&buffer[(char*)&(model->data.options.PsiGS0)-(char*)model]=(FLANDMARK_PSIG*)&buffer[address_offset];
            break;
            case 1: PsiGi = model->data.options.PsiGS1;
                *(FLANDMARK_PSIG**)&buffer[(char*)&(model->data.options.PsiGS1)-(char*)model]=(FLANDMARK_PSIG*)&buffer[address_offset];
            break;
            case 2: PsiGi = model->data.options.PsiGS2;
                *(FLANDMARK_PSIG**)&buffer[(char*)&(model->data.options.PsiGS2)-(char*)model]=(FLANDMARK_PSIG*)&buffer[address_offset];
            break;
        }

        tsize = model->data.options.PSIG_COLS[i]*model->data.options.PSIG_ROWS[i];

        // static part
        disp_pointer = dynamic_offset;
        for (int idx=0; idx < tsize; ++idx)
        {
            address_offset = dynamic_offset;
            memcpy(&buffer[address_offset], &PsiGi[idx], sizeof(FLANDMARK_PSIG));
            dynamic_offset += sizeof(FLANDMARK_PSIG);
        }

        // dynamic part
        for (int idx=0; idx < tsize; ++idx)
        {
            address_offset = dynamic_offset;
            memcpy(&buffer[address_offset], PsiGi[idx].disp, PsiGi[idx].ROWS*PsiGi[idx].COLS*sizeof(int));
            dynamic_offset += PsiGi[idx].ROWS*PsiGi[idx].COLS*sizeof(int);
            // update pointer
            *(int**)&buffer[disp_pointer+idx*sizeof(FLANDMARK_PSIG)+(long long)((char*)&(PsiGi[idx].disp))-(long long)((char*)&(PsiGi[idx]))]=(int*)&buffer[address_offset];
        }
    }

    /* CHECK MODEL  */

    mexPrintf("Checking model...");

    FLANDMARK_Model *tmp = (FLANDMARK_Model*)buffer;

    EError_T err =  flandmark_check_model(model, tmp);

    if (err == NO_ERR)
    {
        mexPrintf(" passed.\n");
    } else {
        mexPrintf(" not passed!\n");

        switch (err)
        {
            case NO_ERR:
                mexErrMsgTxt("No error.\n");
                break;
            case ERROR_M:
                mexErrMsgTxt("model->M does not match!\n");
                break;
            case ERROR_BW:
                mexErrMsgTxt("model->data.bw does not match!\n");
                break;
            case ERROR_BW_MARGIN:
                mexErrMsgTxt("model->data.bw_margin does not match!\n");
                break;
            case ERROR_W:
                mexErrMsgTxt("model->W does not match!\n");
                break;
            case ERROR_DATA_IMAGES:
                mexErrMsgTxt("model->data.Images does not match!\n");
                break;
            case ERROR_DATA_MAPTABLE:
                mexErrMsgTxt("model->data.mapTable does not match!\n");
                break;
            case ERROR_DATA_LBP:
                mexErrMsgTxt("model->data.lbp does not match!\n");
                break;
            case ERROR_DATA_OPTIONS_S:
                mexErrMsgTxt("model->data.options.S does not match!\n");
                break;
            case ERROR_DATA_OPTIONS_PSIG:
                mexErrMsgTxt("model->data.options.PsiGSi does not match!\n");
                break;
            case UNKNOWN_ERROR:
                mexErrMsgTxt("Uknown error!\n");
                break;
        }
    }

    // clean model
    mexPrintf("Cleaning structure model... ");
    flandmark_free(model);
    mexPrintf(" done.\n");
}
