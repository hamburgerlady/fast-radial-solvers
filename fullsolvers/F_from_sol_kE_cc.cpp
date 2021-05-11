#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if (nrhs != 5) {
		mexErrMsgIdAndTxt("F_from_sol_kE_cc:nrhs", "Wrong number of inputs.");
	}
	if (nlhs != 1) {
		mexErrMsgIdAndTxt("F_from_sol_kE_cc:nlhs", "One output required.");
	}

  double *Ar = mxGetPr(prhs[0]);
  double *Br = mxGetPr(prhs[1]);
  double *kk = mxGetPr(prhs[2]);
  double *ff1 = mxGetPr(prhs[3]);
  double *ff2 = mxGetPr(prhs[4]);
  size_t nsol = mxGetM(prhs[2]);
  mwSize i;
  double k,f1,f2;
  
  plhs[0] = mxCreateDoubleMatrix(9,(mwSize)nsol,mxREAL);    
  double *F = mxGetPr(plhs[0]);
  
  for (i=0;i<nsol;i++){
      k = kk[i];
      f1 = ff1[i];
      f2 = ff2[i];
      F[0 + 9*i] = Ar[12]+Br[12]*k+f2*(Ar[0]+Br[0]*k)+f1*(Ar[6]+Br[6]*k);
      F[1 + 9*i] = Ar[13]+Br[13]*k+f2*(Ar[1]+Br[1]*k)+f1*(Ar[7]+Br[7]*k);
      F[2 + 9*i] = f2;
      F[3 + 9*i] = Ar[14]+Br[14]*k+f2*(Ar[2]+Br[2]*k)+f1*(Ar[8]+Br[8]*k);
      F[4 + 9*i] = Ar[15]+Br[15]*k+f2*(Ar[3]+Br[3]*k)+f1*(Ar[9]+Br[9]*k);
      F[5 + 9*i] = f1;
      F[6 + 9*i] = Ar[16]+Br[16]*k+f2*(Ar[4]+Br[4]*k)+f1*(Ar[10]+Br[10]*k);
      F[7 + 9*i] = Ar[17]+Br[17]*k+f2*(Ar[5]+Br[5]*k)+f1*(Ar[11]+Br[11]*k);
      F[8 + 9*i] = 1.0;
  }
}
