#include "mex.h"
#include "math.h"

#define DEG  16
#include "sturm_mart.h"




void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if (nrhs != 1) {
		mexErrMsgIdAndTxt("automatic_generator_cvpr:new_kEf:nrhs", "One input required.");
	}
	if (nlhs != 1) {
		mexErrMsgIdAndTxt("automatic_generator_cvpr:new_kEf:nlhs", "One output required.");
	}    
	if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0])) {
		mexErrMsgIdAndTxt("automatic_generator_cvpr:new_kEf:notDouble", "Input data must be type double.");
	}
	if(mxGetNumberOfElements(prhs[0]) % 17 != 0) {
		mexErrMsgIdAndTxt("automatic_generator_cvpr:new_kEf:incorrectSize", "Input size must be multiple of 17.");
	}
	
    double *p = mxGetPr(prhs[0]);
	plhs[0] = mxCreateDoubleMatrix(1,16,mxREAL);
    double *roots = mxGetPr(plhs[0]);
	int nroots;
    nroots = realRoots (p, roots);
	
	
}


