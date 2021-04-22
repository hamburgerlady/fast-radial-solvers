#include <Eigen/Dense>
#include "mex.h"

#define DEG  32
#include "sturm_mart.h"
#include "charpoly.h"



using namespace Eigen;

void fast_eigenvector_solver(double * eigv, int neig, Eigen::Matrix<double,32,32> &AM, Matrix<std::complex<double>,2,32> &sols);



MatrixXcd solver_new_kEf(const VectorXd& data)
{
	// Compute coefficients
    const double* d = data.data();
    VectorXd coeffs(149);
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
    coeffs[26] = d[32];
    coeffs[27] = d[26];
    coeffs[28] = d[33];
    coeffs[29] = d[40];
    coeffs[30] = d[27];
    coeffs[31] = d[34];
    coeffs[32] = d[41];
    coeffs[33] = d[48];
    coeffs[34] = d[28];
    coeffs[35] = d[35];
    coeffs[36] = d[42];
    coeffs[37] = d[49];
    coeffs[38] = d[56];
    coeffs[39] = d[29];
    coeffs[40] = d[36];
    coeffs[41] = d[43];
    coeffs[42] = d[50];
    coeffs[43] = d[57];
    coeffs[44] = d[30];
    coeffs[45] = d[37];
    coeffs[46] = d[44];
    coeffs[47] = d[51];
    coeffs[48] = d[58];
    coeffs[49] = d[31];
    coeffs[50] = d[38];
    coeffs[51] = d[45];
    coeffs[52] = d[52];
    coeffs[53] = d[59];
    coeffs[54] = d[39];
    coeffs[55] = d[46];
    coeffs[56] = d[53];
    coeffs[57] = d[60];
    coeffs[58] = d[47];
    coeffs[59] = d[54];
    coeffs[60] = d[61];
    coeffs[61] = d[55];
    coeffs[62] = d[62];
    coeffs[63] = d[63];
    coeffs[64] = d[64];
    coeffs[65] = d[65];
    coeffs[66] = d[72];
    coeffs[67] = d[66];
    coeffs[68] = d[73];
    coeffs[69] = d[80];
    coeffs[70] = d[67];
    coeffs[71] = d[74];
    coeffs[72] = d[81];
    coeffs[73] = d[88];
    coeffs[74] = d[68];
    coeffs[75] = d[75];
    coeffs[76] = d[82];
    coeffs[77] = d[89];
    coeffs[78] = d[96];
    coeffs[79] = d[69];
    coeffs[80] = d[76];
    coeffs[81] = d[83];
    coeffs[82] = d[90];
    coeffs[83] = d[97];
    coeffs[84] = d[70];
    coeffs[85] = d[77];
    coeffs[86] = d[84];
    coeffs[87] = d[91];
    coeffs[88] = d[98];
    coeffs[89] = d[71];
    coeffs[90] = d[78];
    coeffs[91] = d[85];
    coeffs[92] = d[92];
    coeffs[93] = d[99];
    coeffs[94] = d[79];
    coeffs[95] = d[86];
    coeffs[96] = d[93];
    coeffs[97] = d[100];
    coeffs[98] = d[87];
    coeffs[99] = d[94];
    coeffs[100] = d[101];
    coeffs[101] = d[95];
    coeffs[102] = d[102];
    coeffs[103] = d[103];
    coeffs[104] = d[104];
    coeffs[105] = d[105];
    coeffs[106] = d[113];
    coeffs[107] = d[106];
    coeffs[108] = d[114];
    coeffs[109] = d[122];
    coeffs[110] = d[107];
    coeffs[111] = d[115];
    coeffs[112] = d[123];
    coeffs[113] = d[131];
    coeffs[114] = d[108];
    coeffs[115] = d[116];
    coeffs[116] = d[124];
    coeffs[117] = d[132];
    coeffs[118] = d[140];
    coeffs[119] = d[109];
    coeffs[120] = d[117];
    coeffs[121] = d[125];
    coeffs[122] = d[133];
    coeffs[123] = d[141];
    coeffs[124] = d[110];
    coeffs[125] = d[118];
    coeffs[126] = d[126];
    coeffs[127] = d[134];
    coeffs[128] = d[142];
    coeffs[129] = d[111];
    coeffs[130] = d[119];
    coeffs[131] = d[127];
    coeffs[132] = d[135];
    coeffs[133] = d[143];
    coeffs[134] = d[112];
    coeffs[135] = d[120];
    coeffs[136] = d[128];
    coeffs[137] = d[136];
    coeffs[138] = d[144];
    coeffs[139] = d[121];
    coeffs[140] = d[129];
    coeffs[141] = d[137];
    coeffs[142] = d[145];
    coeffs[143] = d[130];
    coeffs[144] = d[138];
    coeffs[145] = d[146];
    coeffs[146] = d[139];
    coeffs[147] = d[147];
    coeffs[148] = d[148];



	// Setup elimination template
	static const int coeffs0_ind[] = { 0,24,64,104,5,2,29,69,109,9,5,33,73,113,13,8,37,9,77,33,5,73,117,9,38,78,118,13,43,83,38,9,78,123,61,21,101,146,23,63,23,103,148,1,25,0,65,24,64,105,2,0,26,66,106,8,4,32,5,72,29,2,69,112,17,12,42,13,82,37,8,77,122,9,5,17,48,88,43,13,83,128,9 };
	static const int coeffs1_ind[] = { 3,27,1,67,25,65,107,0,4,1,28,2,68,26,0,66,108,6,30,3,70,27,67,110,0,1,7,3,31,4,71,28,1,68,111,2,0,10,34,6,74,30,70,114,1,3,11,6,35,7,75,31,0,3,71,115,2,4,1,12,7,36,8,76,32,4,72,116,5,2,14,39,10,79,34,74,119,3,6,15,10,40,11,80,35,1,6,75,120,4,7,3,16,11,41,12,81,36,2,7,76,121,5,8,4,44,14,84,39,79,124,6,10,18,14,45,15,85,40,3,10,80,125,7,11,6,19,15,46,16,86,41,4,11,81,126,8,12,7,20,16,47,17,87,42,5,12,82,127,9,13,8,49,89,44,84,129,10,14,50,18,90,45,6,14,85,130,11,15,10,21,18,51,19,91,46,7,15,86,131,12,16,11,22,19,52,20,92,47,8,16,87,132,13,17,12,20,53,93,48,9,17,88,133,13,49,89,134,14,54,94,50,10,90,135,15,18,14,55,21,95,51,11,18,91,136,16,19,15,23,21,56,22,96,52,12,19,92,137,17,20,16,22,57,97,53,13,20,93,138,17,54,14,94,139,18,58,98,55,15,95,140,19,21,18,59,23,99,56,16,21,96,141,20,22,19,23,60,100,57,17,22,97,142,20,58,18,98,143,21,61,101,59,19,99,144,22,23,21,62,102,60,20,23,100,145,22,63,103,62,22,102,147,23 };
	static const int C0_ind[] = { 0,2,4,9,13,14,15,17,22,26,27,28,30,35,39,40,41,42,43,44,46,47,48,53,54,56,61,66,67,69,70,72,73,74,83,84,86,87,88,96,97,99,100,104,106,107,108,109,112,113,117,118,119,121,126,130,131,132,133,134,135,137,138,139,143,144,145,146,147,148,150,151,152,154,155,157,158,160,161,163,164,165,168 } ;
	static const int C1_ind[] = { 0,2,3,4,5,8,9,11,13,14,15,16,17,18,20,21,22,26,28,29,30,31,34,35,36,37,39,40,41,42,43,44,46,47,48,50,51,52,54,55,56,57,60,61,62,63,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,85,86,87,89,90,91,93,94,95,96,99,100,101,102,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,132,133,134,135,138,139,140,141,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,184,186,187,190,191,192,193,197,198,199,200,201,202,203,204,205,206,207,208,209,210,211,212,213,214,215,216,217,218,219,220,221,222,223,224,225,226,227,228,229,230,231,232,233,235,236,238,239,240,241,242,243,246,252,255,256,257,262,264,265,266,268,269,270,271,272,275,276,277,278,279,280,281,282,283,284,285,286,287,288,289,290,291,292,293,294,295,296,297,298,300,301,303,304,305,306,307,308,311,317,318,320,321,322,327,329,330,331,333,334,335,336,337,340,341,342,343,344,345,346,347,348,349,350,352,353,355,356,357,358,359,360,363,369,370,372,373,374,379,381,382,383,385,386,387,388,389,392,394,395,396,397,398,399,402,405,407,408,409,411,412,415 };

	Matrix<double,13,13> C0; C0.setZero();
	Matrix<double,13,32> C1; C1.setZero();
	for (int i = 0; i < 83; i++) { C0(C0_ind[i]) = coeffs(coeffs0_ind[i]); }
	for (int i = 0; i < 314; i++) { C1(C1_ind[i]) = coeffs(coeffs1_ind[i]); } 

	Matrix<double,13,32> C12 = C0.partialPivLu().solve(C1);




	// Setup action matrix
	Matrix<double,37, 32> RR;
	RR << -C12.bottomRows(5), Matrix<double,32,32>::Identity(32, 32);

	static const int AM_ind[] = { 0,1,5,6,7,8,2,9,10,11,12,13,14,3,15,16,17,18,4,19,20,21,22,23,25,26,27,28,30,31,32,35 };
	Matrix<double, 32, 32> AM;
	for (int i = 0; i < 32; i++) {
		AM.row(i) = RR.row(AM_ind[i]);
	}

	Matrix<std::complex<double>, 2, 32> sols;
	sols.setZero();

	// Solve eigenvalue problem


	double p[1+32];
	Matrix<double, 32, 32> AMp = AM;
	charpoly_danilevsky_piv(AMp, p);	
	double roots[32];
	int nroots;
	// find_real_roots_sturm(p, 32, roots, &nroots, 8, 0);
    nroots = realRoots (p, roots);
	fast_eigenvector_solver(roots, nroots, AM, sols);






	return sols;
}
// Action =  y
// Quotient ring basis (V) = x^4*y^6,x^3*y^7,x^4*y^5,x^3*y^6,x^4*y^4,x^3*y^5,x^2*y^6,x^4*y^3,x^3*y^4,x^2*y^5,x^4*y^2,x^3*y^3,x^2*y^4,x*y^5,x^4*y,x^3*y^2,x^2*y^3,x*y^4,y^5,x^4,x^3*y,x^2*y^2,x*y^3,y^4,x^3,x^2*y,x*y^2,y^3,x^2,x*y,y^2,y,
// Available monomials (RR*V) = x^4*y^7,x^3*y^8,x^2*y^7,x*y^6,y^6,x^4*y^6,x^3*y^7,x^4*y^5,x^3*y^6,x^4*y^4,x^3*y^5,x^2*y^6,x^4*y^3,x^3*y^4,x^2*y^5,x^4*y^2,x^3*y^3,x^2*y^4,x*y^5,x^4*y,x^3*y^2,x^2*y^3,x*y^4,y^5,x^4,x^3*y,x^2*y^2,x*y^3,y^4,x^3,x^2*y,x*y^2,y^3,x^2,x*y,y^2,y,


void fast_eigenvector_solver(double * eigv, int neig, Eigen::Matrix<double,32,32> &AM, Matrix<std::complex<double>,2,32> &sols) {
	static const int ind[] = { 0,1,6,13,18 };	
	// Truncated action matrix containing non-trivial rows
	Matrix<double, 5, 32> AMs;
	double zi[8];
	
	for (int i = 0; i < 5; i++)	{
		AMs.row(i) = AM.row(ind[i]);
	}
	for (int i = 0; i < neig; i++) {
		zi[0] = eigv[i];
		for (int j = 1; j < 8; j++)
		{
			zi[j] = zi[j - 1] * eigv[i];
		}
		Matrix<double, 5,5> AA;
        AA.col(0) = zi[5] * AMs.col(0) + zi[4] * AMs.col(2) + zi[3] * AMs.col(4) + zi[2] * AMs.col(7) + zi[1] * AMs.col(10) + zi[0] * AMs.col(14) + AMs.col(19);
        AA.col(1) = zi[6] * AMs.col(1) + zi[5] * AMs.col(3) + zi[4] * AMs.col(5) + zi[3] * AMs.col(8) + zi[2] * AMs.col(11) + zi[1] * AMs.col(15) + zi[0] * AMs.col(20) + AMs.col(24);
        AA.col(2) = zi[5] * AMs.col(6) + zi[4] * AMs.col(9) + zi[3] * AMs.col(12) + zi[2] * AMs.col(16) + zi[1] * AMs.col(21) + zi[0] * AMs.col(25) + AMs.col(28);
        AA.col(3) = zi[4] * AMs.col(13) + zi[3] * AMs.col(17) + zi[2] * AMs.col(22) + zi[1] * AMs.col(26) + zi[0] * AMs.col(29);
        AA.col(4) = zi[4] * AMs.col(18) + zi[3] * AMs.col(23) + zi[2] * AMs.col(27) + zi[1] * AMs.col(30) + zi[0] * AMs.col(31);
        AA(0,0) = AA(0,0) - zi[6];
        AA(1,1) = AA(1,1) - zi[7];
        AA(2,2) = AA(2,2) - zi[6];
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
		mexErrMsgIdAndTxt("automatic_generator_cvpr:new_kEf:nrhs", "One input required.");
	}
	if (nlhs != 1) {
		mexErrMsgIdAndTxt("automatic_generator_cvpr:new_kEf:nlhs", "One output required.");
	}    
	if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0])) {
		mexErrMsgIdAndTxt("automatic_generator_cvpr:new_kEf:notDouble", "Input data must be type double.");
	}
	if(mxGetNumberOfElements(prhs[0]) % 149 != 0) {
		mexErrMsgIdAndTxt("automatic_generator_cvpr:new_kEf:incorrectSize", "Input size must be multiple of 149.");
	}
	int n_instances = mxGetNumberOfElements(prhs[0]) / 149;
	double *input = mxGetPr(prhs[0]);
	plhs[0] = mxCreateDoubleMatrix(2,32*n_instances,mxCOMPLEX);
	double* zr = mxGetPr(plhs[0]);
	double* zi = mxGetPi(plhs[0]);
	for(int k = 0; k < n_instances; k++) {
		const VectorXd data = Map<const VectorXd>(input + k*149, 149);
		MatrixXcd sols = solver_new_kEf(data);
		Index offset = k*sols.size();
		for (Index i = 0; i < sols.size(); i++) {
			zr[i+offset] = sols(i).real();
			zi[i+offset] = sols(i).imag();
		}
	}
}


