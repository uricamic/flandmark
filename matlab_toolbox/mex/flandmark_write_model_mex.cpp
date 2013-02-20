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
    if (nrhs < 2 || nrhs > 3)
    {
        mexErrMsgTxt("Improper number of input arguments.\n\n"
                "FLANDMARK_WRITE_MODEL loads model for flandmark detector and writes it to file.\n\n"
                "Synopsis: \n"
                "flandmark_write_model(model, fname)\n"
                "\n"
                "Input: \n"
                "  model [structure]\n"
                "  fname [string]\n"
            );
    }

    // variables for interface to matlab
    const mwSize *dims, *tmp_dims;
    mxArray *field_ptr, *field_data_ptr, *tmp_field_ptr, *cell_ptr;
    // pointers to data types used in mat file
    double * p_double = 0;
    uint32_t * p_uint32 = 0;
    int tsize = -1, tmp_tsize = -1;
    EError_T retval = UNKNOWN_ERROR;
    char * fname = mxArrayToString(prhs[1]);

    // check input
    if (!mxIsStruct(prhs[0]))
    {
        mexErrMsgTxt("Improper argument type \"model\".\n\n"
                "  model must be a structure containing joint parameter vector W and structure data.\n"
             );
    }
    mxArray *in_w = mxGetField(prhs[0], 0, "W");
    mxArray *in_data = mxGetField(prhs[0], 0, "data");
    if (in_w == NULL || in_data == NULL)
    {
        mexErrMsgTxt("Improper argument type \"model\".\n\n"
                "  model must be a structure containing joint parameter vector \"W\" and structure \"data\".\n"
             );
    }
    if (fname == NULL)
    {
        mexErrMsgTxt("Improper argument type \"fname\"\n\n"
                "  fname must be a string.\n"
              );
    }


    // Model structure (see detector_model.h for details)
    FLANDMARK_Model * model = (FLANDMARK_Model*)malloc(sizeof(FLANDMARK_Model));
    uint8_t M = 0;
    // model.data.options.M --------------------------------------------------
    mexPrintf("Getting number of components M:\n");
    field_data_ptr = mxGetField(prhs[0], 0, "data");
    field_ptr = mxGetField(field_data_ptr, 0, "options");
    tmp_field_ptr = mxGetField(field_ptr, 0, "M");
    p_double = (double*)mxGetPr(tmp_field_ptr);
    model->data.options.M = (uint8_t)p_double[0];
    M = model->data.options.M;
    mexPrintf("model->data.options.M = %d\n", model->data.options.M);
    // model->data.options.bw -------------------------------------------------
    tmp_field_ptr = mxGetField(field_ptr, 0, "bw");
    p_double = (double*)mxGetPr(tmp_field_ptr);
    model->data.options.bw[0] = (int)p_double[0];
    model->data.options.bw[1] = (int)p_double[1];
    mexPrintf("model->data.options.bw = [%d, %d]\n",
            model->data.options.bw[0], model->data.options.bw[1]);
    model->data.imSize[0] = model->data.options.bw[0];
    model->data.imSize[1] = model->data.options.bw[1];
    // model->data.options.bw_margin ------------------------------------------
    tmp_field_ptr = mxGetField(field_ptr, 0, "bw_margin");
    p_double = (double*)mxGetPr(tmp_field_ptr);
    model->data.options.bw_margin[0] = (int)p_double[0];
    model->data.options.bw_margin[1] = (int)p_double[1];
    mexPrintf("model->data.options.bw_margin = [%d, %d]\n",
            model->data.options.bw_margin[0], model->data.options.bw_margin[1]);
    // model->W ---------------------------------------------------------------
    // get model->W from mat file
    mexPrintf("Loading field model->W:\n");
    field_ptr = mxGetField(prhs[0], 0, "W");
    dims = mxGetDimensions(field_ptr);
    model->W_ROWS = dims[0]; model->W_COLS = dims[1];
    mexPrintf("model->W = [%d x %d]\n", model->W_ROWS, model->W_COLS);
    p_double = (double*)mxGetPr(field_ptr);
    model->W = (double*)malloc(model->W_ROWS * sizeof(double));
    for (int i = 0; i < model->W_ROWS; ++i)
    {
        model->W[i] = p_double[i];
    }
    // model->data.mapTable ----------------------------------------------------
    mexPrintf("Loading field model->data.mapTable:\n");
    field_ptr = mxGetField(field_data_ptr, 0, "mapTable");
    dims = mxGetDimensions(field_ptr);
    mexPrintf("mapTable = [%d x %d]\n", dims[0], dims[1]);
    mexPrintf(" = %d\n", dims[0]*dims[1]);
    if (dims[0]*dims[1] != model->data.options.M*4)
    {
        mexErrMsgTxt("Error. Possibly corrupted file.\n");
    }
    model->data.mapTable = (int*)calloc(M*4, sizeof(int));
    p_double = (double*)mxGetPr(field_ptr);
    mexPrintf("model->data.mapTable: \n");
    for (int i = 0; i < model->data.options.M*4; ++i)
    {
        model->data.mapTable[i] = (int)p_double[i];
        mexPrintf(" %d", model->data.mapTable[i]);
    }
    mexPrintf("\n");

    // model->data.lbp ---------------------------------------------------------
    mexPrintf("Loading field model->data.lbp:\n");
    field_ptr = mxGetField(field_data_ptr, 0, "lbp");
    dims = mxGetDimensions(field_ptr);
    mexPrintf("lbp = [%d x %d]\n", dims[0], dims[1]);
    model->data.lbp = (FLANDMARK_LBP*)malloc(M*sizeof(FLANDMARK_LBP));
    for (int idx = 0; idx < model->data.options.M; ++idx)
    {
        // get pointer to cell lbp{idx}
        cell_ptr = mxGetCell(field_ptr, idx);
        // get lbp{idx}.winSize
        tmp_field_ptr = mxGetField(cell_ptr, 0, "winSize");
        tmp_dims = mxGetDimensions(tmp_field_ptr);
        mexPrintf("lbp[%d].winSize = [%d x %d]", idx, tmp_dims[0], tmp_dims[1]);
        p_double = (double*)mxGetPr(tmp_field_ptr);
        for (int i = 0; i < 2; ++i)
        {
            model->data.lbp[idx].winSize[i] = (int)p_double[i];
            mexPrintf(" %d", model->data.lbp[idx].winSize[i]);
        }
        mexPrintf("\n");

        // get lbp{idx}.hop
        model->data.lbp[idx].hop = 4;

        // get lbp{idx}.wins
        tmp_field_ptr = mxGetField(cell_ptr, 0, "wins");
        tmp_dims = mxGetDimensions(tmp_field_ptr);
        model->data.lbp[idx].WINS_ROWS = tmp_dims[0]; model->data.lbp[idx].WINS_COLS = tmp_dims[1];
        mexPrintf("lbp[%d].wins = [%d x %d]\n", idx, model->data.lbp[idx].WINS_ROWS, model->data.lbp[idx].WINS_COLS);
        p_uint32 = (uint32_t*)mxGetPr(tmp_field_ptr);
        tsize = model->data.lbp[idx].WINS_ROWS * model->data.lbp[idx].WINS_COLS;
        model->data.lbp[idx].wins = (uint32_t*)malloc(tsize * sizeof(uint32_t));
        for (int i = 0; i < tsize; ++i)
        {
            model->data.lbp[idx].wins[i] = p_uint32[i];
        }
    }

    // model->data.options.S ---------------------------------------------------
    mexPrintf("Loading field model->data.options.S:\n");
    //field_data_ptr = mxGetField(prhs[0], 0, "data");
    field_ptr = mxGetField(field_data_ptr, 0, "options");
    tmp_field_ptr = mxGetField(field_ptr, 0, "S");
    dims = mxGetDimensions(tmp_field_ptr);
    mexPrintf("options.S = [%d x %d]\n", dims[0], dims[1]);
    model->data.options.S = (int*)calloc(4*M, sizeof(int));
    if (dims[0]*dims[1] != 4*model->data.options.M)
    {
        mexErrMsgTxt("Error. Possibly corrupted file.\n");
    }
    p_double = (double*)mxGetPr(tmp_field_ptr);
    for (int i = 0; i < 4*model->data.options.M; ++i)
    {
        model->data.options.S[i] = (int)p_double[i];
        mexPrintf(" %d", model->data.options.S[i]);
    }
    mexPrintf("\n");

    // model->data.options.PsiG ------------------------------------------------
    mexPrintf("Loading field model->data.options.PsiG:\n");
    field_ptr = mxGetField(field_data_ptr, 0, "options");

    FLANDMARK_PSIG * PsiGi = NULL;
    char psig_id[7] = "PsiGS#";

    for (int psig_idx = 0; psig_idx < 3; ++psig_idx)
    {
        mexPrintf("PsiGS%d:\n", psig_idx);

        switch (psig_idx)
        {
            case 0:
                //PsiGi = model->data.options.PsiGS0;
                tmp_field_ptr = mxGetField(field_ptr, 0, "PsiGS0");
                psig_id[5] = '0';
                break;
            case 1:
                //PsiGi = model->data.options.PsiGS1;
                tmp_field_ptr = mxGetField(field_ptr, 0, "PsiGS1");
                psig_id[5] = '1';
                break;
            case 2:
                //PsiGi = model->data.options.PsiGS2;
                tmp_field_ptr = mxGetField(field_ptr, 0, "PsiGS2");
                psig_id[5] = '2';
                break;
        }

        if (tmp_field_ptr == NULL)
        {
            mexErrMsgTxt("Looks like you are trying to load a model using gtab version of deformation cost.\n"
                    "We now support only loading models with gdisp deformation cost.\n"
                    "\n"
                 );
        }
        dims = mxGetDimensions(tmp_field_ptr);
        tsize = dims[0]*dims[1];
        model->data.options.PSIG_ROWS[psig_idx] = dims[0]; model->data.options.PSIG_COLS[psig_idx] = dims[1];

        mexPrintf("options.%s = [%d x %d]\n", psig_id, model->data.options.PSIG_ROWS[psig_idx], model->data.options.PSIG_COLS[psig_idx]);

        switch (psig_idx)
        {
            case 0:
                model->data.options.PsiGS0 = (FLANDMARK_PSIG*)malloc(tsize*sizeof(FLANDMARK_PSIG));
                PsiGi = model->data.options.PsiGS0;
                break;
            case 1:
                model->data.options.PsiGS1 = (FLANDMARK_PSIG*)malloc(tsize*sizeof(FLANDMARK_PSIG));
                PsiGi = model->data.options.PsiGS1;
                break;
            case 2:
                model->data.options.PsiGS2 = (FLANDMARK_PSIG*)malloc(tsize*sizeof(FLANDMARK_PSIG));
                PsiGi = model->data.options.PsiGS2;
                break;
        }

        for (int idx = 0; idx < tsize; ++idx)
        {
            cell_ptr = mxGetCell(tmp_field_ptr, idx);
            tmp_dims = mxGetDimensions(cell_ptr);
            tmp_tsize = tmp_dims[0]*tmp_dims[1];
            PsiGi[idx].ROWS = tmp_dims[0]; PsiGi[idx].COLS = tmp_dims[1];
            PsiGi[idx].disp = (int*)mxGetPr(cell_ptr);
        }
    }

    model->bb = 0;
    model->sf = 0;
    model->normalizedImageFrame = 0;

    mexPrintf("Structure model loaded.\n");

    // write model
    mexPrintf("\n\n");
    mexPrintf("Calling flandmark_writeModel...\n");
    flandmark_write_model(fname, model);
    mexPrintf("Write model done.\n");

    // read model
    mexPrintf("Calling flandmark_init...\n");
    FLANDMARK_Model * tst = flandmark_init(fname);
    mexPrintf("Read model done.\n");

    // chceck model
    mexPrintf("Calling flandmark_checkModel... ");
    retval = flandmark_check_model(model, tst);
    switch (retval)
    {
        case NO_ERR:
            mexPrintf("all done.\n");
        break;
        case ERROR_M:
            mexErrMsgTxt("check model->data.options.M not passed.\n");
        break;
        case ERROR_BW:
            mexErrMsgTxt("check model->data.options.bw not passed.\n");
        break;
        case ERROR_BW_MARGIN:
            mexErrMsgTxt("check model->data.options.bw_margin not passed.\n");
        break;

        case ERROR_W:
            mexErrMsgTxt("check model->w not passed.\n");
        break;
        case ERROR_DATA_IMAGES:
            mexErrMsgTxt("check model->data.Images not passed.\n");
        break;
        case ERROR_DATA_MAPTABLE:
            mexErrMsgTxt("check model->data.mapTable not passed.\n");
        break;

        case ERROR_DATA_LBP:
            mexErrMsgTxt("check model->data.lbp not passed.\n");
        break;
        case ERROR_DATA_OPTIONS_S:
            mexErrMsgTxt("check model->data.options.S not passed.\n");
        break;
        case ERROR_DATA_OPTIONS_PSIG:
            mexErrMsgTxt("check model->data.options.PsiG not passed.\n");
        break;

        default:
            mexErrMsgTxt("Unknown error.\n");
        break;
    }

    // cleanup
    mexPrintf("Cleanup...");
    flandmark_free(model);
    flandmark_free(tst);
    mexPrintf(" done.\n");
}
