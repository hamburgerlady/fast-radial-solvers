#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if (nrhs != 4) {
		mexErrMsgIdAndTxt("F_from_sol_kfEf_cc:nrhs", "Wrong number of inputs.");
	}
	if (nlhs != 1) {
		mexErrMsgIdAndTxt("F_from_sol_kfEf_cc:nlhs", "One output required.");
	}

  double *Ar = mxGetPr(prhs[0]);
  double *Br = mxGetPr(prhs[1]);
  double *kk = mxGetPr(prhs[2]);
  double *ff1 = mxGetPr(prhs[3]);
  size_t nsol = mxGetM(prhs[2]);
  mwSize i;
  double k,f1;
  
  plhs[0] = mxCreateDoubleMatrix(9,(mwSize)nsol,mxREAL);    
  double *F = mxGetPr(plhs[0]);
  
  for (i=0;i<nsol;i++){
      k = kk[i];
      f1 = ff1[i];
 
      F[0 + 9*i] = Ar[0]*Ar[20]-Ar[6]*Ar[14]+Ar[0]*Ar[13]*f1-Ar[6]*Ar[7]*f1+Ar[0]*Br[20]*k-Ar[6]*Br[14]*k-Ar[14]*Br[6]*k+Ar[20]*Br[0]*k+Br[0]*Br[20]*(k*k)-Br[6]*Br[14]*(k*k)+Br[0]*Br[13]*f1*(k*k)-Br[6]*Br[7]*f1*(k*k)+Ar[0]*Br[13]*f1*k-Ar[6]*Br[7]*f1*k-Ar[7]*Br[6]*f1*k+Ar[13]*Br[0]*f1*k;
      F[1 + 9*i] = Ar[1]*Ar[20]-Ar[6]*Ar[15]+Ar[1]*Ar[13]*f1-Ar[6]*Ar[8]*f1+Ar[1]*Br[20]*k-Ar[6]*Br[15]*k-Ar[15]*Br[6]*k+Ar[20]*Br[1]*k+Br[1]*Br[20]*(k*k)-Br[6]*Br[15]*(k*k)+Br[1]*Br[13]*f1*(k*k)-Br[6]*Br[8]*f1*(k*k)+Ar[1]*Br[13]*f1*k-Ar[6]*Br[8]*f1*k-Ar[8]*Br[6]*f1*k+Ar[13]*Br[1]*f1*k;
      F[2 + 9*i] = -Ar[20]-Ar[13]*f1-Br[20]*k-Br[13]*f1*k;
      F[3 + 9*i] = Ar[2]*Ar[20]-Ar[6]*Ar[16]+Ar[2]*Ar[13]*f1-Ar[6]*Ar[9]*f1+Ar[2]*Br[20]*k-Ar[6]*Br[16]*k-Ar[16]*Br[6]*k+Ar[20]*Br[2]*k+Br[2]*Br[20]*(k*k)-Br[6]*Br[16]*(k*k)+Br[2]*Br[13]*f1*(k*k)-Br[6]*Br[9]*f1*(k*k)+Ar[2]*Br[13]*f1*k-Ar[6]*Br[9]*f1*k-Ar[9]*Br[6]*f1*k+Ar[13]*Br[2]*f1*k;
      F[4 + 9*i] = Ar[3]*Ar[20]-Ar[6]*Ar[17]+Ar[3]*Ar[13]*f1-Ar[6]*Ar[10]*f1+Ar[3]*Br[20]*k-Ar[6]*Br[17]*k-Ar[17]*Br[6]*k+Ar[20]*Br[3]*k+Br[3]*Br[20]*(k*k)-Br[6]*Br[17]*(k*k)+Br[3]*Br[13]*f1*(k*k)-Br[6]*Br[10]*f1*(k*k)+Ar[3]*Br[13]*f1*k-Ar[6]*Br[10]*f1*k-Ar[10]*Br[6]*f1*k+Ar[13]*Br[3]*f1*k;
      F[5 + 9*i] = f1*(Ar[6]+Br[6]*k);
      F[6 + 9*i] = Ar[4]*Ar[20]-Ar[6]*Ar[18]+Ar[4]*Ar[13]*f1-Ar[6]*Ar[11]*f1+Ar[4]*Br[20]*k-Ar[6]*Br[18]*k-Ar[18]*Br[6]*k+Ar[20]*Br[4]*k+Br[4]*Br[20]*(k*k)-Br[6]*Br[18]*(k*k)+Br[4]*Br[13]*f1*(k*k)-Br[6]*Br[11]*f1*(k*k)+Ar[4]*Br[13]*f1*k-Ar[6]*Br[11]*f1*k-Ar[11]*Br[6]*f1*k+Ar[13]*Br[4]*f1*k;
      F[7 + 9*i] = Ar[5]*Ar[20]-Ar[6]*Ar[19]+Ar[5]*Ar[13]*f1-Ar[6]*Ar[12]*f1+Ar[5]*Br[20]*k-Ar[6]*Br[19]*k-Ar[19]*Br[6]*k+Ar[20]*Br[5]*k+Br[5]*Br[20]*(k*k)-Br[6]*Br[19]*(k*k)+Br[5]*Br[13]*f1*(k*k)-Br[6]*Br[12]*f1*(k*k)+Ar[5]*Br[13]*f1*k-Ar[6]*Br[12]*f1*k-Ar[12]*Br[6]*f1*k+Ar[13]*Br[5]*f1*k;
      F[8 + 9*i] = Ar[6]+Br[6]*k;
          
      
  }
}
