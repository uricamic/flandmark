#include "mex.h" // 使用MEX文件必须包含的头文件
// 执行具体工作的C函数
double add(double x, double y)
{
    return x + y;
}
// MEX文件接口函数
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