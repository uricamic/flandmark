#include "mex.h" // ʹ��MEX�ļ����������ͷ�ļ�
// ִ�о��幤����C����
double add(double x, double y)
{
    return x + y;
}
// MEX�ļ��ӿں���
void mexFunction(int nlhs,mxArray *plhs[], int nrhs,const mxArray *prhs[])
{
    double *a;
    double b, c;
    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    a = mxGetPr(plhs[0]);
    b = *(mxGetPr(prhs[0]));
    c = *(mxGetPr(prhs[1]));
    *a = add(b, c);
}