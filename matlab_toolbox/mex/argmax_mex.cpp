#include "flandmark_detector.h"

#include <mex.h>
#include <cstring>
#include <cfloat>

void get_psi_sparse(FLANDMARK_PSI_SPARSE* Psi, const mxArray *matlab_psi, int idx)
{
    mxArray *cell_ptr = mxGetCell(matlab_psi, idx);
    const mwSize *dims;

    dims = mxGetDimensions(cell_ptr);

    Psi->PSI_COLS = dims[1];
    Psi->PSI_ROWS = dims[0];
    Psi->idxs = (uint32_t*)mxGetPr(cell_ptr);
}

void mexFunction( int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[] )
{
    if(nrhs < 4 || nrhs > 5)
     mexErrMsgTxt("Improper number of input arguments.\n\n"
                  "ARGMAX_MEX estimates the position of landmarks.\n\n"
                  "Synopsis: \n"
                  "  [ \\hat{s} ]= argmax_mex(options, W, mapTable, lbp [, L]) \n"
                  "\n"
                  "  Input: \n"
                  "    options [1 x 1 structure] detector options\n"
                  "    W [n x 1 (double)] joint parameter vector\n"
                  "    mapTable [M x 4 (double)] mapping table for joint parameter vector W\n"
                  "    lbp [1 x M (cell)] LBP sparse features\n"
                  "    L [M x 1 (cell)] precomputed losses, optional parameter (for learning stage)\n"
                  "  Output: \n"
                  "    \\hat{s} [M x 2 (double)] estimated landmark positions\n");
    
    //--------------------------------------------------------------------------
    // Load input
    double * W;
    uint32_t W_ROWS, W_COLS;
    
    FLANDMARK_Options options;
    uint8_t M;
    const mwSize *dims, *tmp_dims;
    double * tmp_ptr, * p_double;
    mxArray *tmp_field_ptr, *cell_ptr;
    int tsize = -1, tmp_tsize = -1, *p_int;
    
    mxArray * psigs0 = mxGetField(prhs[0], 0, "PsiGS0");
    dims = mxGetDimensions(psigs0);
    options.PSIG_ROWS[0] = dims[0];
    options.PSIG_COLS[0] = dims[1];
    mxArray * psigs1 = mxGetField(prhs[0], 0, "PsiGS1");
    dims = mxGetDimensions(psigs1);
    options.PSIG_ROWS[1] = dims[0];
    options.PSIG_COLS[1] = dims[1];
    mxArray * psigs2 = mxGetField(prhs[0], 0, "PsiGS2");
    dims = mxGetDimensions(psigs2);
    options.PSIG_ROWS[2] = dims[0];
    options.PSIG_COLS[2] = dims[1];

    // model->data.options.PsiG ------------------------------------------------
    FLANDMARK_PSIG * PsiGi = NULL;

    for (int psig_idx = 0; psig_idx < 3; ++psig_idx)
    {
        switch (psig_idx)
        {
            case 0:
                tmp_field_ptr = psigs0;
                break;
            case 1:
                tmp_field_ptr = psigs1;
                break;
            case 2:
                tmp_field_ptr = psigs2;
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
        options.PSIG_ROWS[psig_idx] = dims[0];
        options.PSIG_COLS[psig_idx] = dims[1];

        switch (psig_idx)
        {
            case 0:
                options.PsiGS0 = (FLANDMARK_PSIG*)malloc(tsize*sizeof(FLANDMARK_PSIG));
                PsiGi = options.PsiGS0;
                break;
            case 1:
                options.PsiGS1 = (FLANDMARK_PSIG*)malloc(tsize*sizeof(FLANDMARK_PSIG));
                PsiGi = options.PsiGS1;
                break;
            case 2:
                options.PsiGS2 = (FLANDMARK_PSIG*)malloc(tsize*sizeof(FLANDMARK_PSIG));
                PsiGi = options.PsiGS2;
                break;
        }

        for (int idx = 0; idx < tsize; ++idx)
        {
            cell_ptr = mxGetCell(tmp_field_ptr, idx);
            tmp_dims = mxGetDimensions(cell_ptr);
            tmp_tsize = tmp_dims[0]*tmp_dims[1];
            PsiGi[idx].ROWS = tmp_dims[0];
            PsiGi[idx].COLS = tmp_dims[1];
            PsiGi[idx].disp = (int*)mxGetPr(cell_ptr);
        }
    }

    tmp_field_ptr = mxGetField(prhs[0], 0, "bw");
    tmp_ptr = (double*)mxGetPr(tmp_field_ptr);
    options.bw[0] = (int)tmp_ptr[0];
    options.bw[1] = (int)tmp_ptr[1];
    tmp_field_ptr = mxGetField(prhs[0], 0, "bw_margin");
    tmp_ptr = (double*)mxGetPr(tmp_field_ptr);
    options.bw_margin[0] = (int)tmp_ptr[0];
    options.bw_margin[1] = (int)tmp_ptr[1];

    int rows = -1, cols = -1;
    double * result;    
    bool lossy = (nrhs > 4) ? true : false;
    
    // options ----------------------------------------------------------------
    //options.M = 7;
    tmp_field_ptr = mxGetField(prhs[0], 0, "M");
    tmp_ptr = (double*)mxGetPr(tmp_field_ptr);
    options.M = (int)tmp_ptr[0];
    M = options.M;
    
    //int mapTable[7*4];
    int * mapTable = (int*)calloc(M*4, sizeof(int));
    options.S = (int*)calloc(4*M, sizeof(int));
    
    // options.S ---------------------------------------------------
    tmp_field_ptr = mxGetField(prhs[0], 0, "S");
    dims = mxGetDimensions(tmp_field_ptr);
    tmp_ptr = (double*)mxGetPr(tmp_field_ptr);
    for (int i = 0; i < 4*options.M; ++i)
    {
        options.S[i] = (int)tmp_ptr[i];
    }
    
    // W ----------------------------------------------------------------------
    dims = mxGetDimensions(prhs[1]);
    W_ROWS = dims[0]; W_COLS = dims[1];
    bool WisSparse = mxIsSparse(prhs[1]);
    if (WisSparse)
    {
        W = (double*)calloc(W_ROWS*W_COLS, sizeof(double));
        double * pr = (double*)mxGetPr(prhs[1]);
        mwIndex * ir = mxGetIr(prhs[1]);
        mwSize nzmax = mxGetNzmax(prhs[1]);
        for (int ii = 0; ii < nzmax; ++ii)
        {
            W[ir[ii]] = pr[ii];
        }
    } else {
        W = (double*)mxGetPr(prhs[1]);
    }
    // ------------------------------------------------------------------------

    // mapTable
    tmp_ptr = (double*)mxGetPr(prhs[2]);
    dims = mxGetDimensions(prhs[2]);
    if (dims[0]*dims[1] != M*4)
    {
        mexErrMsgTxt("mapTable must be a matrix [M x 4 (double)]");
    }
    for (int i = 0; i < M*4; ++i)
    {
        mapTable[i] = (int)tmp_ptr[i];
    }
    
    // prepare output matrix
    plhs[0] = mxCreateNumericMatrix(2, M, mxDOUBLE_CLASS, mxREAL);
    result = (double*)mxGetPr(plhs[0]);
    
    // get PSI matrix
    FLANDMARK_PSI_SPARSE * Psi_sparse = (FLANDMARK_PSI_SPARSE*)malloc(M*sizeof(FLANDMARK_PSI_SPARSE));
    for (int idx = 0; idx < M; ++idx)
    {
        get_psi_sparse(&Psi_sparse[idx], prhs[3], idx);
    }

    //--------------------------------------------------------------------------
    // get Q, G
    double ** q = (double**)calloc(M, sizeof(double*));
    double ** g = (double**)calloc((M-1), sizeof(double*));

    int idx_qtemp = 0;

    double * L;
    
    for (int idx = 0; idx < M; ++idx)
    {
        // Q
        tsize = mapTable[INDEX(idx, 1, M)] - mapTable[INDEX(idx, 0, M)] + 1;

        double * q_temp = (double*)calloc(tsize, sizeof(double));
        memcpy(q_temp, W+mapTable[INDEX(idx, 0, M)]-1, tsize*sizeof(double));

        // sparse dot product <W_q, PSI_q>
        cell_ptr = mxGetCell(prhs[3], idx);
        dims = mxGetDimensions(cell_ptr);
        cols = dims[1];
        rows = dims[0];
        uint32_t *psi_temp = (uint32_t*)mxGetPr(cell_ptr);
        if (lossy)
        {
            L = (double*)mxGetPr(mxGetCell(prhs[4], idx));
        }
        q[idx] = (double*)malloc(cols*sizeof(double));
        for (int i = 0; i < cols; ++i)
        {
            double dotprod = 0.0f;
            for (int j = 0; j < rows; ++j)
            {
                idx_qtemp = psi_temp[(rows*i) + j];
                dotprod += q_temp[ idx_qtemp ];
            }
            q[idx][i] = dotprod;
            if (lossy)
            {
                q[idx][i] += L[i];
            }
        }
        free(q_temp);

        // G
        if (idx > 0)
        {
            tsize = mapTable[INDEX(idx, 3, M)] - mapTable[INDEX(idx, 2, M)] + 1;
            g[idx - 1] = (double*)malloc(tsize*sizeof(double));
            memcpy(g[idx - 1], W+mapTable[INDEX(idx, 2, M)]-1, tsize*sizeof(double));
        }
    }

    // compute argmax
    flandmark_argmax(result, &options, mapTable, Psi_sparse, q, g);

    // cleanup W if it was created
    if (WisSparse)
    {
        free(W);
    }
    
    // cleanup q
    for (int i = 0; i < M; ++i)
    {
        free(q[i]);
    }
    free(q);

    // cleanup g
    for (int i = 0; i < M - 1; ++i)
    {
        free(g[i]);
    }
    free(g);

    free(options.S);
    free(options.PsiGS0);
    free(options.PsiGS1);
    free(options.PsiGS2);
    free(mapTable);
    free(Psi_sparse);
    
    return;
}
