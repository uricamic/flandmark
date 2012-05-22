#ifndef __FLANDMARK_DETECTOR_H_
#define __FLANDMARK_DETECTOR_H_

#include <stdint.h>
#include <cv.h>
#include <cvaux.h>

// index row-order matrices
#define INDEX(ROW, COL, NUM_ROWS) ((COL)*(NUM_ROWS)+(ROW))
#define ROW(IDX, ROWS) (((IDX)-1) % (ROWS))
#define COL(IDX, ROWS) (((IDX)-1) / (ROWS))

typedef struct psig_struct {
    int * disp;
    int ROWS, COLS;
} FLANDMARK_PSIG;

typedef struct options_struct {
    uint8_t M;
    int * S;
    int bw[2], bw_margin[2];
    FLANDMARK_PSIG *PsiGS0, *PsiGS1, *PsiGS2;
    int PSIG_ROWS[3], PSIG_COLS[3];
} FLANDMARK_Options;

typedef struct lbp_struct {
    int winSize[2];
    uint8_t hop;
    uint32_t *wins;
    int WINS_ROWS, WINS_COLS;
} FLANDMARK_LBP;

typedef struct data_struct {
    FLANDMARK_LBP * lbp;
    int Images_ROWS, Images_COLS;
    int imSize[2];
    int * mapTable;
    FLANDMARK_Options options;
} FLANDMARK_Data;

typedef struct model_struct {
    double * W;
    int W_ROWS, W_COLS;
    FLANDMARK_Data data;
    uint8_t *normalizedImageFrame;
    IplImage *croppedImage;
    IplImage *resizedImage;
    double *bb;
    float *sf;
    int p_width, p_height;
} FLANDMARK_Model;

typedef struct psi_struct {
    char * data;
    uint32_t PSI_ROWS, PSI_COLS;
} FLANDMARK_PSI;

typedef struct psi_sparse {
    uint32_t * idxs;
    uint32_t PSI_ROWS, PSI_COLS;
} FLANDMARK_PSI_SPARSE;
// -------------------------------------------------------------------------

enum EError_T {
    NO_ERR=0,
    ERROR_M=1,
    ERROR_BW=2,
    ERROR_BW_MARGIN=3,
    ERROR_W=4,
    ERROR_DATA_IMAGES=5,
    ERROR_DATA_MAPTABLE=6,
    ERROR_DATA_LBP=7,
    ERROR_DATA_OPTIONS_S=8,
    ERROR_DATA_OPTIONS_PSIG=9,
    UNKNOWN_ERROR=100,
};

// read / write structure Model from / to file procedures

/**
 * Function flandmark_init
 *
 * Given the path to the file containing the model in binary form, this function will return a pointer to this model. It returns null pointer in the case of failure
 *
 * \param[in] filename
 * \return Pointer to the FLANDMARK_Model data structure
 */
FLANDMARK_Model * flandmark_init(const char* filename);

/**
 * Function flandmark_write model
 *
 * This function writes given FLANDMARK_model data structure to a file specified by its path.
 *
 * \param[in] filename
 * \param[in] model
 */
void flandmark_writeModel(const char* filename, FLANDMARK_Model* model);

/**
 * Function flandmark_checkModel
 *
 * This function checks if both given FLANDMARK_Model data structres are equal
 *
 * \param[in] model
 * \param[in] tst
 * \return
 */
EError_T flandmark_checkModel(FLANDMARK_Model* model, FLANDMARK_Model* tst);

/**
 * Function flandmark_free
 *
 * This function dealocates the FLANDMARK_Model data structure
 *
 * \param[in] model
 */
void flandmark_free(FLANDMARK_Model* model);

// getPsiMat (calls LBP features computation - liblbpfeatures from LIBOCAS)
/**
 *
 * \param[out] Psi
 * \param[in] model
 * \param[in] lbpidx
 */
void flandmark_getPsiMat(FLANDMARK_PSI* Psi, FLANDMARK_Model* model, int lbpidx);

/**
 * Computes LBP features representing it as sparse matrix (i.e. only inices with ones are stored in connected list)
 *
 * \param[out] Psi
 * \param[in] model
 * \param[in] lbpidx
 */
void flandmark_getPsiMatSparse(FLANDMARK_PSI_SPARSE* Psi, FLANDMARK_Model* model, int lbpidx);

// dot product maximization with max-index return
/**
 * Function maximizedotprod
 *
 * \param[in]
 * \param[in]
 * \param[out]
 */
void maximize_gdotprod(double *maximum, double *idx, double *first, double *second, int *third, int cols, int tsize);

/**
 * Function getNormalizedImageFrame
 *
 *
 */
int getNormalizedImageFrame(IplImage *input, const int bbox[], double *bb, uint8_t *face_img, FLANDMARK_Model *model);

int imCrop(IplImage *input, IplImage *output, const CvRect region, FLANDMARK_Model *model);

int imResize(IplImage *input, IplImage *output, int width, int height, FLANDMARK_Model *model);

/**
 * Function flandmark_detect_base
 *
 * Estimates positions of facial landmarks in the normalized image frame.
 *
 * \param[in] face_image pointer to 1D uint8 array with normalized image frame of face
 * \param[in] model Data structure holding info about model
 * \param[in, out] int array representing 2D array of size [2 x options.M] with estimated positions of landmarks
 * \return int indicator of success or fail of the detection
 */
int flandmark_detect_base(uint8_t *face_image, FLANDMARK_Model *model, float *landmarks);

/**
 * Function flandmark_detect
 *
 * Estimates positions of facial landmarks given the image and the bounding box of the detected face
 *
 */
int flandmark_detect(IplImage *img, int * bbox, FLANDMARK_Model *model, float * landmarks, int * bw_margin = 0);

#endif // __LIBFLD_DETECTOR_H_
