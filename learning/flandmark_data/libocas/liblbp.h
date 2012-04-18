#ifndef _liblbp_h
#define _liblbp_h

#include <stdint.h>

#define LIBLBP_INDEX(ROW,COL,NUM_ROWS) ((COL)*(NUM_ROWS)+(ROW))
#define LIBLBP_MIN(A,B) ((A) > (B) ? (B) : (A))

extern void liblbp_pyr_features(char *vec, uint32_t vec_nDim, uint32_t *img, uint16_t img_nRows, uint16_t img_nCols );
extern double liblbp_pyr_dotprod(double *vec, uint32_t vec_nDim, uint32_t *img, uint16_t img_nRows, uint16_t img_nCols);
extern void liblbp_pyr_addvec(int64_t *vec, uint32_t vec_nDim, uint32_t *img, uint16_t img_nRows, uint16_t img_nCols);
extern void liblbp_pyr_subvec(int64_t *vec, uint32_t vec_nDim, uint32_t *img, uint16_t img_nRows, uint16_t img_nCols);
extern uint32_t liblbp_pyr_get_dim(uint16_t img_nRows, uint16_t img_nCols, uint16_t nPyramids);

#endif
