#include "argmax_mex.h"

#include <mex.h>
#include <cstring>
#include <cfloat>

void maximize_gdotprod(double * maximum, double * idx, double * first, double * second, double * third, int cols, int tsize)
{
    *maximum = -FLT_MAX;
    *idx = -1;
    double *temp_psig = (double*)calloc(tsize, sizeof(double));
    for (int dp_i = 0; dp_i < cols; ++dp_i)
    {
        memcpy(temp_psig, third+(dp_i*tsize), tsize*sizeof(double));
        double dotprod = 0.0f;
        for (int dp_j = 0; dp_j < tsize; ++dp_j)
        {
            dotprod += second[dp_j]*temp_psig[dp_j];
        }
        if (*maximum < first[dp_i]+dotprod)
        {
            *idx = dp_i;
            *maximum = first[dp_i]+dotprod;
        }
    }
    free(temp_psig);
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
    const mwSize *dims;
    double * tmp_ptr;
    mxArray *tmp_field_ptr, *cell_ptr;
    
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

    int tsize = -1, rows = -1, cols = -1;
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
    int * indices = (int*)malloc(M*sizeof(int));
    tsize = mapTable[INDEX(1, 3, M)] - mapTable[INDEX(1, 2, M)] + 1;

    // left branch - store maximum and index of s5 for all positions of s1
    dims = mxGetDimensions(mxGetCell(prhs[3], 1)); 
    int q1_length = dims[1];
    
    double * s1 = (double *)calloc(2*q1_length, sizeof(double));
    double * s1_maxs = (double *)calloc(q1_length, sizeof(double));
    for (int i = 0; i < q1_length; ++i)
    {
        // dot product <g_5, PsiGS1>
        dims = mxGetDimensions(mxGetCell(psigs1, INDEX(i, 0, options.PSIG_ROWS[1])));
        maximize_gdotprod(
            &s1[INDEX(0, i, 2)], &s1[INDEX(1, i, 2)],
            q[5], g[4], (double*)mxGetPr(mxGetCell(psigs1, INDEX(i, 0, options.PSIG_ROWS[1]))),
            dims[1], tsize);
        
        s1[INDEX(0, i, 2)] += q[1][i];
    }
    for (int i = 0; i < q1_length; ++i)
    {
        s1_maxs[i] = s1[INDEX(0, i, 2)];
    }
    
    // right branch (s2->s6) - store maximum and index of s6 for all positions of s2
    dims = mxGetDimensions(mxGetCell(prhs[3], 2));
    int q2_length = dims[1];
    
    double * s2 = (double *)calloc(2*q2_length, sizeof(double));
    double * s2_maxs = (double *)calloc(q2_length, sizeof(double));
    for (int i = 0; i < q2_length; ++i)
    {
        // dot product <g_6, PsiGS2>     
        dims = mxGetDimensions(mxGetCell(psigs2, INDEX(i, 0, options.PSIG_ROWS[2])));
        maximize_gdotprod(
                &s2[INDEX(0, i, 2)], (double*)&s2[INDEX(1, i, 2)],
                q[6], g[5], (double*)mxGetPr(mxGetCell(psigs2, INDEX(i, 0, options.PSIG_ROWS[2]))),
                dims[1], tsize);
        s2[INDEX(0, i, 2)] += q[2][i];
    }
    for (int i = 0; i < q2_length; ++i)
    {
        s2_maxs[i] = s2[INDEX(0, i, 2)];
    }
   
    // the root s0 and its connections
    dims = mxGetDimensions(mxGetCell(prhs[3], 0)); 
    int q0_length = dims[1];
    double maxs0 = -FLT_MAX; 
    int maxs0_idx = -1;
    double maxq10 = -FLT_MAX, maxq20 = -FLT_MAX, maxq30 = -FLT_MAX, maxq40 = -FLT_MAX, 
            maxq70 = -FLT_MAX;
    double * s0 = (double *)calloc(M*q0_length, sizeof(double));
    
    for (int i = 0; i < q0_length; ++i)
    {
        // q10
        maxq10 = -FLT_MAX;
        dims = mxGetDimensions(mxGetCell(psigs0, INDEX(i, 0, options.PSIG_ROWS[0])));
        maximize_gdotprod(
                &maxq10, &s0[INDEX(1, i, M)],
                s1_maxs, g[0], (double*)mxGetPr(mxGetCell(psigs0, INDEX(i, 0, options.PSIG_ROWS[0]))),
                dims[1], tsize);
        s0[INDEX(5, i, M)] = s1[INDEX(1, (int)s0[INDEX(1, i, M)], 2)];
        
        // q20
        maxq20 = -FLT_MAX;
        dims = mxGetDimensions(mxGetCell(psigs0, INDEX(i, 1, options.PSIG_ROWS[0])));
        maximize_gdotprod(
                &maxq20, &s0[INDEX(2, i, M)],
                s2_maxs, g[1], (double*)mxGetPr(mxGetCell(psigs0, INDEX(i, 1, options.PSIG_ROWS[0]))),
                dims[1], tsize);
        s0[INDEX(6, i, M)] = s2[INDEX(1, (int)s0[INDEX(2, i, M)], 2)];
        
        // q30
        maxq30 = -FLT_MAX;
        dims = mxGetDimensions(mxGetCell(psigs0, INDEX(i, 2, options.PSIG_ROWS[0])));
        maximize_gdotprod(
                &maxq30, &s0[INDEX(3, i, M)],
                q[3], g[2], (double*)mxGetPr(mxGetCell(psigs0, INDEX(i, 2, options.PSIG_ROWS[0]))),
                dims[1], tsize);
        
        // q40
        maxq40 = -FLT_MAX;
        dims = mxGetDimensions(mxGetCell(psigs0, INDEX(i, 3, options.PSIG_ROWS[0])));
        maximize_gdotprod(
                &maxq40, &s0[INDEX(4, i, M)],
                q[4], g[3], (double*)mxGetPr(mxGetCell(psigs0, INDEX(i, 3, options.PSIG_ROWS[0]))),
                dims[1], tsize);
        
        // q70
        maxq70 = -FLT_MAX;
        dims = mxGetDimensions(mxGetCell(psigs0, INDEX(i, 4, options.PSIG_ROWS[0])));
        maximize_gdotprod(
                &maxq70, &s0[INDEX(7, i, M)],
                q[7], g[6], (double*)mxGetPr(mxGetCell(psigs0, INDEX(i, 4, options.PSIG_ROWS[0]))),
                dims[1], tsize);
        
        // sum q10+q20+q30+q40+q70
        if (maxs0 < maxq10+maxq20+maxq30+maxq40+maxq70+q[0][i])
        {
            maxs0_idx = i;
            s0[INDEX(0, i, M)] = i;
            maxs0 = maxq10+maxq20+maxq30+maxq40+maxq70+q[0][i];
        }
    }
    
    // get indices
    for (int i = 0; i < M; ++i)
    {
        indices[i] = (int)s0[INDEX(0, maxs0_idx, M)+i]+1;
    }

    // cleanup temp variables
    free(s0);
    free(s1); free(s1_maxs);
    free(s2); free(s2_maxs);

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
    
    // convert 1D indices to 2D coordinates of estimated positions
    int * optionsS = &options.S[0];
    for (int i = 0; i < M; ++i)
    {
        int rows = optionsS[INDEX(3, i, 4)] - optionsS[INDEX(1, i, 4)] + 1;
        result[INDEX(0, i, 2)] = COL(indices[i], rows) + optionsS[INDEX(0, i, 4)];
        result[INDEX(1, i, 2)] = ROW(indices[i], rows) + optionsS[INDEX(1, i, 4)];
    }
    free(indices);
    free(options.S);
    free(mapTable);
    
    return;
}
