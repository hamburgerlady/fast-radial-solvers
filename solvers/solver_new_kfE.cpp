#include <Eigen/Dense>
#include "mex.h"

#define DEG  27
#include "sturm_mart.h"
#include "charpoly.h"



using namespace Eigen;

void fast_eigenvector_solver(double * eigv, int neig, Eigen::Matrix<double,27,27> &AM, Matrix<std::complex<double>,2,27> &sols);



MatrixXcd solver_new_kfE(const VectorXd& data)
{
	// Compute coefficients
    const double* d = data.data();
    VectorXd coeffs(124);
    coeffs[0] = d[0];
    coeffs[1] = d[1];
    coeffs[2] = d[6];
    coeffs[3] = d[2];
    coeffs[4] = d[7];
    coeffs[5] = d[12];
    coeffs[6] = d[3];
    coeffs[7] = d[8];
    coeffs[8] = d[13];
    coeffs[9] = d[18];
    coeffs[10] = d[4];
    coeffs[11] = d[9];
    coeffs[12] = d[14];
    coeffs[13] = d[19];
    coeffs[14] = d[5];
    coeffs[15] = d[10];
    coeffs[16] = d[15];
    coeffs[17] = d[20];
    coeffs[18] = d[11];
    coeffs[19] = d[16];
    coeffs[20] = d[21];
    coeffs[21] = d[17];
    coeffs[22] = d[22];
    coeffs[23] = d[23];
    coeffs[24] = d[24];
    coeffs[25] = d[25];
    coeffs[26] = d[30];
    coeffs[27] = d[26];
    coeffs[28] = d[31];
    coeffs[29] = d[36];
    coeffs[30] = d[27];
    coeffs[31] = d[32];
    coeffs[32] = d[37];
    coeffs[33] = d[42];
    coeffs[34] = d[28];
    coeffs[35] = d[33];
    coeffs[36] = d[38];
    coeffs[37] = d[43];
    coeffs[38] = d[48];
    coeffs[39] = d[29];
    coeffs[40] = d[34];
    coeffs[41] = d[39];
    coeffs[42] = d[44];
    coeffs[43] = d[49];
    coeffs[44] = d[35];
    coeffs[45] = d[40];
    coeffs[46] = d[45];
    coeffs[47] = d[50];
    coeffs[48] = d[41];
    coeffs[49] = d[46];
    coeffs[50] = d[51];
    coeffs[51] = d[47];
    coeffs[52] = d[52];
    coeffs[53] = d[53];
    coeffs[54] = d[54];
    coeffs[55] = d[55];
    coeffs[56] = d[61];
    coeffs[57] = d[56];
    coeffs[58] = d[62];
    coeffs[59] = d[68];
    coeffs[60] = d[57];
    coeffs[61] = d[63];
    coeffs[62] = d[69];
    coeffs[63] = d[75];
    coeffs[64] = d[58];
    coeffs[65] = d[64];
    coeffs[66] = d[70];
    coeffs[67] = d[76];
    coeffs[68] = d[82];
    coeffs[69] = d[59];
    coeffs[70] = d[65];
    coeffs[71] = d[71];
    coeffs[72] = d[77];
    coeffs[73] = d[83];
    coeffs[74] = d[60];
    coeffs[75] = d[66];
    coeffs[76] = d[72];
    coeffs[77] = d[78];
    coeffs[78] = d[84];
    coeffs[79] = d[67];
    coeffs[80] = d[73];
    coeffs[81] = d[79];
    coeffs[82] = d[85];
    coeffs[83] = d[74];
    coeffs[84] = d[80];
    coeffs[85] = d[86];
    coeffs[86] = d[81];
    coeffs[87] = d[87];
    coeffs[88] = d[88];
    coeffs[89] = d[89];
    coeffs[90] = d[90];
    coeffs[91] = d[96];
    coeffs[92] = d[91];
    coeffs[93] = d[97];
    coeffs[94] = d[103];
    coeffs[95] = d[92];
    coeffs[96] = d[98];
    coeffs[97] = d[104];
    coeffs[98] = d[110];
    coeffs[99] = d[93];
    coeffs[100] = d[99];
    coeffs[101] = d[105];
    coeffs[102] = d[111];
    coeffs[103] = d[117];
    coeffs[104] = d[94];
    coeffs[105] = d[100];
    coeffs[106] = d[106];
    coeffs[107] = d[112];
    coeffs[108] = d[118];
    coeffs[109] = d[95];
    coeffs[110] = d[101];
    coeffs[111] = d[107];
    coeffs[112] = d[113];
    coeffs[113] = d[119];
    coeffs[114] = d[102];
    coeffs[115] = d[108];
    coeffs[116] = d[114];
    coeffs[117] = d[120];
    coeffs[118] = d[109];
    coeffs[119] = d[115];
    coeffs[120] = d[121];
    coeffs[121] = d[116];
    coeffs[122] = d[122];
    coeffs[123] = d[123];



	// Setup elimination template
	static const int coeffs0_ind[] = { 0,2,26,56,91,1,0,4,2,26,28,58,93,2,5,29,59,94,3,1,7,4,28,31,61,96,4,2,8,5,29,32,62,97,0,24,54,89,5,9,33,63,98,9,38,68,103 };
	static const int coeffs1_ind[] = { 1,0,24,25,55,90,3,1,25,27,57,92,6,3,27,30,60,95,10,6,30,34,64,99,6,3,11,7,31,35,65,100,7,4,12,8,32,36,66,101,8,5,13,9,33,37,67,102,14,10,34,39,69,104,10,6,15,11,35,40,70,105,11,7,16,12,36,41,71,106,12,8,17,13,37,42,72,107,13,9,38,43,73,108,14,39,74,109,14,10,18,15,40,44,75,110,15,11,19,16,41,45,76,111,16,12,20,17,42,46,77,112,17,13,43,47,78,113,14,18,44,79,114,18,15,21,19,45,48,80,115,19,16,22,20,46,49,81,116,20,17,47,50,82,117,18,21,48,83,118,21,19,23,22,49,51,84,119,22,20,50,52,85,120,21,23,51,86,121,23,22,52,53,87,122,23,53,88,123 };
	static const int C0_ind[] = { 0,2,5,6,7,8,9,10,11,12,13,14,15,16,18,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,42,45,46,47,48,50,53,54,55,56,61,62,63 } ;
	static const int C1_ind[] = { 2,3,4,5,6,7,10,11,12,13,14,15,18,19,20,21,22,23,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,92,93,94,95,99,100,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,132,133,134,135,137,139,140,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,164,165,166,167,169,171,172,174,175,176,177,178,179,180,181,182,183,184,185,188,189,190,191,193,195,196,198,199,200,201,204,205,206,207,209,212,214,215 };

	Matrix<double,8,8> C0; C0.setZero();
	Matrix<double,8,27> C1; C1.setZero();
	for (int i = 0; i < 47; i++) { C0(C0_ind[i]) = coeffs(coeffs0_ind[i]); }
	for (int i = 0; i < 179; i++) { C1(C1_ind[i]) = coeffs(coeffs1_ind[i]); } 

	Matrix<double,8,27> C12 = C0.partialPivLu().solve(C1);




	// Setup action matrix
	Matrix<double,32, 27> RR;
	RR << -C12.bottomRows(5), Matrix<double,27,27>::Identity(27, 27);

	static const int AM_ind[] = { 2,5,6,7,0,1,3,8,9,10,11,4,12,13,14,15,16,18,19,20,21,23,24,25,27,28,30 };
	Matrix<double, 27, 27> AM;
	for (int i = 0; i < 27; i++) {
		AM.row(i) = RR.row(AM_ind[i]);
	}

	Matrix<std::complex<double>, 2, 27> sols;
	sols.setZero();

	// Solve eigenvalue problem


	double p[1+27];
	Matrix<double, 27, 27> AMp = AM;
	charpoly_danilevsky_piv(AMp, p);	
	double roots[27];
	int nroots;
	// find_real_roots_sturm(p, 27, roots, &nroots, 8, 0);
    nroots = realRoots (p, roots);
	fast_eigenvector_solver(roots, nroots, AM, sols);






	return sols;
}
// Action =  y
// Quotient ring basis (V) = x^4*y^5,x^4*y^4,x^4*y^3,x^4*y^2,x^3*y^3,x^2*y^4,x*y^5,x^4*y,x^3*y^2,x^2*y^3,x*y^4,y^5,x^4,x^3*y,x^2*y^2,x*y^3,y^4,x^3,x^2*y,x*y^2,y^3,x^2,x*y,y^2,x,y,1,
// Available monomials (RR*V) = x^3*y^4,x^2*y^5,x^4*y^6,x*y^6,y^6,x^4*y^5,x^4*y^4,x^4*y^3,x^4*y^2,x^3*y^3,x^2*y^4,x*y^5,x^4*y,x^3*y^2,x^2*y^3,x*y^4,y^5,x^4,x^3*y,x^2*y^2,x*y^3,y^4,x^3,x^2*y,x*y^2,y^3,x^2,x*y,y^2,x,y,1,


void fast_eigenvector_solver(double * eigv, int neig, Eigen::Matrix<double,27,27> &AM, Matrix<std::complex<double>,2,27> &sols) {
	static const int ind[] = { 0,4,5,6,11 };	
	// Truncated action matrix containing non-trivial rows
	Matrix<double, 5, 27> AMs;
	double zi[6];
	
	for (int i = 0; i < 5; i++)	{
		AMs.row(i) = AM.row(ind[i]);
	}
	for (int i = 0; i < neig; i++) {
		zi[0] = eigv[i];
		for (int j = 1; j < 6; j++)
		{
			zi[j] = zi[j - 1] * eigv[i];
		}
		Matrix<double, 5,5> AA;
        AA.col(0) = zi[4] * AMs.col(0) + zi[3] * AMs.col(1) + zi[2] * AMs.col(2) + zi[1] * AMs.col(3) + zi[0] * AMs.col(7) + AMs.col(12);
        AA.col(1) = zi[2] * AMs.col(4) + zi[1] * AMs.col(8) + zi[0] * AMs.col(13) + AMs.col(17);
        AA.col(2) = zi[3] * AMs.col(5) + zi[2] * AMs.col(9) + zi[1] * AMs.col(14) + zi[0] * AMs.col(18) + AMs.col(21);
        AA.col(3) = zi[4] * AMs.col(6) + zi[3] * AMs.col(10) + zi[2] * AMs.col(15) + zi[1] * AMs.col(19) + zi[0] * AMs.col(22) + AMs.col(24);
        AA.col(4) = zi[4] * AMs.col(11) + zi[3] * AMs.col(16) + zi[2] * AMs.col(20) + zi[1] * AMs.col(23) + zi[0] * AMs.col(25) + AMs.col(26);
        AA(0,0) = AA(0,0) - zi[5];
        AA(1,1) = AA(1,1) - zi[3];
        AA(2,2) = AA(2,2) - zi[4];
        AA(3,3) = AA(3,3) - zi[5];
        AA(4,4) = AA(4,4) - zi[5];


		Matrix<double, 4, 1>  s = AA.leftCols(4).colPivHouseholderQr().solve(-AA.col(4));
        sols(1,i) = s(3);
        sols(0,i) = zi[0];

	}
}



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if (nrhs != 1) {
		mexErrMsgIdAndTxt("automatic_generator_cvpr:new_kfE:nrhs", "One input required.");
	}
	if (nlhs != 1) {
		mexErrMsgIdAndTxt("automatic_generator_cvpr:new_kfE:nlhs", "One output required.");
	}    
	if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0])) {
		mexErrMsgIdAndTxt("automatic_generator_cvpr:new_kfE:notDouble", "Input data must be type double.");
	}
	if(mxGetNumberOfElements(prhs[0]) % 124 != 0) {
		mexErrMsgIdAndTxt("automatic_generator_cvpr:new_kfE:incorrectSize", "Input size must be multiple of 124.");
	}
	int n_instances = mxGetNumberOfElements(prhs[0]) / 124;
	double *input = mxGetPr(prhs[0]);
	plhs[0] = mxCreateDoubleMatrix(2,27*n_instances,mxCOMPLEX);
	double* zr = mxGetPr(plhs[0]);
	double* zi = mxGetPi(plhs[0]);
	for(int k = 0; k < n_instances; k++) {
		const VectorXd data = Map<const VectorXd>(input + k*124, 124);
		MatrixXcd sols = solver_new_kfE(data);
		Index offset = k*sols.size();
		for (Index i = 0; i < sols.size(); i++) {
			zr[i+offset] = sols(i).real();
			zi[i+offset] = sols(i).imag();
		}
	}
}


