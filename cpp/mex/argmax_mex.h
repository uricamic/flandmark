#ifndef __ARGMAX_MEX_H__
#define __ARGMAX_MEX_H__

#include <stdint.h>

#define INDEX(ROW,COL,NUM_ROWS) ((COL)*(NUM_ROWS)+(ROW))
#define ROW(IDX, ROWS) (((IDX)-1) % (ROWS))
#define COL(IDX, ROWS) (((IDX)-1) / (ROWS))

typedef struct psi_sparse {
    uint32_t * idxs;
    uint32_t PSI_ROWS, PSI_COLS;
} FLANDMARK_PSI_SPARSE;

typedef struct options_struct {
    uint8_t M;
    //int S[4*7];
    int * S;
    int PSIG_ROWS[3], PSIG_COLS[3];
} FLANDMARK_Options;

// dot product maximization with max-index return
/**
 * Function maximizedotprod
 *
 * \param[in]
 * \param[in]
 * \param[out]
 */
//void maximize_gdotprod(double *maximum, double *idx, double *first, double *second, double *third, int cols, int tsize);
void maximize_gdotprod(double *maximum, double *idx, double *first, double *second, int *third, int cols, int tsize);

#endif //__ARGMAX_MEX_H__