#include <Eigen/Dense>
#include "mex.h"

#define DEG  26
#include "sturm_mart.h"
#include "charpoly.h"



using namespace Eigen;

void fast_eigenvector_solver(double * eigv, int neig, Eigen::Matrix<double,26,26> &AM, Matrix<std::complex<double>,3,26> &sols);



MatrixXcd solver_new_kE(const VectorXd& data)
{
	// Compute coefficients
    const double* d = data.data();
    VectorXd coeffs(360);
    coeffs[0] = d[0];
    coeffs[1] = d[4];
    coeffs[2] = d[12];
    coeffs[3] = d[24];
    coeffs[4] = d[1];
    coeffs[5] = d[5];
    coeffs[6] = d[13];
    coeffs[7] = d[25];
    coeffs[8] = d[8];
    coeffs[9] = d[16];
    coeffs[10] = d[28];
    coeffs[11] = d[2];
    coeffs[12] = d[6];
    coeffs[13] = d[14];
    coeffs[14] = d[26];
    coeffs[15] = d[9];
    coeffs[16] = d[17];
    coeffs[17] = d[29];
    coeffs[18] = d[20];
    coeffs[19] = d[32];
    coeffs[20] = d[3];
    coeffs[21] = d[7];
    coeffs[22] = d[15];
    coeffs[23] = d[27];
    coeffs[24] = d[10];
    coeffs[25] = d[18];
    coeffs[26] = d[30];
    coeffs[27] = d[21];
    coeffs[28] = d[33];
    coeffs[29] = d[36];
    coeffs[30] = d[11];
    coeffs[31] = d[19];
    coeffs[32] = d[31];
    coeffs[33] = d[22];
    coeffs[34] = d[34];
    coeffs[35] = d[37];
    coeffs[36] = d[23];
    coeffs[37] = d[35];
    coeffs[38] = d[38];
    coeffs[39] = d[39];
    coeffs[40] = d[40];
    coeffs[41] = d[44];
    coeffs[42] = d[52];
    coeffs[43] = d[64];
    coeffs[44] = d[41];
    coeffs[45] = d[45];
    coeffs[46] = d[53];
    coeffs[47] = d[65];
    coeffs[48] = d[48];
    coeffs[49] = d[56];
    coeffs[50] = d[68];
    coeffs[51] = d[42];
    coeffs[52] = d[46];
    coeffs[53] = d[54];
    coeffs[54] = d[66];
    coeffs[55] = d[49];
    coeffs[56] = d[57];
    coeffs[57] = d[69];
    coeffs[58] = d[60];
    coeffs[59] = d[72];
    coeffs[60] = d[43];
    coeffs[61] = d[47];
    coeffs[62] = d[55];
    coeffs[63] = d[67];
    coeffs[64] = d[50];
    coeffs[65] = d[58];
    coeffs[66] = d[70];
    coeffs[67] = d[61];
    coeffs[68] = d[73];
    coeffs[69] = d[76];
    coeffs[70] = d[51];
    coeffs[71] = d[59];
    coeffs[72] = d[71];
    coeffs[73] = d[62];
    coeffs[74] = d[74];
    coeffs[75] = d[77];
    coeffs[76] = d[63];
    coeffs[77] = d[75];
    coeffs[78] = d[78];
    coeffs[79] = d[79];
    coeffs[80] = d[80];
    coeffs[81] = d[83];
    coeffs[82] = d[89];
    coeffs[83] = d[98];
    coeffs[84] = d[81];
    coeffs[85] = d[84];
    coeffs[86] = d[90];
    coeffs[87] = d[99];
    coeffs[88] = d[86];
    coeffs[89] = d[92];
    coeffs[90] = d[101];
    coeffs[91] = d[82];
    coeffs[92] = d[85];
    coeffs[93] = d[91];
    coeffs[94] = d[100];
    coeffs[95] = d[87];
    coeffs[96] = d[93];
    coeffs[97] = d[102];
    coeffs[98] = d[95];
    coeffs[99] = d[104];
    coeffs[100] = d[88];
    coeffs[101] = d[94];
    coeffs[102] = d[103];
    coeffs[103] = d[96];
    coeffs[104] = d[105];
    coeffs[105] = d[107];
    coeffs[106] = d[97];
    coeffs[107] = d[106];
    coeffs[108] = d[108];
    coeffs[109] = d[109];
    coeffs[110] = d[110];
    coeffs[111] = d[114];
    coeffs[112] = d[122];
    coeffs[113] = d[134];
    coeffs[114] = d[111];
    coeffs[115] = d[115];
    coeffs[116] = d[123];
    coeffs[117] = d[135];
    coeffs[118] = d[118];
    coeffs[119] = d[126];
    coeffs[120] = d[138];
    coeffs[121] = d[112];
    coeffs[122] = d[116];
    coeffs[123] = d[124];
    coeffs[124] = d[136];
    coeffs[125] = d[119];
    coeffs[126] = d[127];
    coeffs[127] = d[139];
    coeffs[128] = d[130];
    coeffs[129] = d[142];
    coeffs[130] = d[113];
    coeffs[131] = d[117];
    coeffs[132] = d[125];
    coeffs[133] = d[137];
    coeffs[134] = d[120];
    coeffs[135] = d[128];
    coeffs[136] = d[140];
    coeffs[137] = d[131];
    coeffs[138] = d[143];
    coeffs[139] = d[146];
    coeffs[140] = d[121];
    coeffs[141] = d[129];
    coeffs[142] = d[141];
    coeffs[143] = d[132];
    coeffs[144] = d[144];
    coeffs[145] = d[147];
    coeffs[146] = d[133];
    coeffs[147] = d[145];
    coeffs[148] = d[148];
    coeffs[149] = d[149];
    coeffs[150] = d[150];
    coeffs[151] = d[154];
    coeffs[152] = d[162];
    coeffs[153] = d[174];
    coeffs[154] = d[151];
    coeffs[155] = d[155];
    coeffs[156] = d[163];
    coeffs[157] = d[175];
    coeffs[158] = d[158];
    coeffs[159] = d[166];
    coeffs[160] = d[178];
    coeffs[161] = d[152];
    coeffs[162] = d[156];
    coeffs[163] = d[164];
    coeffs[164] = d[176];
    coeffs[165] = d[159];
    coeffs[166] = d[167];
    coeffs[167] = d[179];
    coeffs[168] = d[170];
    coeffs[169] = d[182];
    coeffs[170] = d[153];
    coeffs[171] = d[157];
    coeffs[172] = d[165];
    coeffs[173] = d[177];
    coeffs[174] = d[160];
    coeffs[175] = d[168];
    coeffs[176] = d[180];
    coeffs[177] = d[171];
    coeffs[178] = d[183];
    coeffs[179] = d[186];
    coeffs[180] = d[161];
    coeffs[181] = d[169];
    coeffs[182] = d[181];
    coeffs[183] = d[172];
    coeffs[184] = d[184];
    coeffs[185] = d[187];
    coeffs[186] = d[173];
    coeffs[187] = d[185];
    coeffs[188] = d[188];
    coeffs[189] = d[189];
    coeffs[190] = d[190];
    coeffs[191] = d[193];
    coeffs[192] = d[199];
    coeffs[193] = d[208];
    coeffs[194] = d[191];
    coeffs[195] = d[194];
    coeffs[196] = d[200];
    coeffs[197] = d[209];
    coeffs[198] = d[196];
    coeffs[199] = d[202];
    coeffs[200] = d[211];
    coeffs[201] = d[192];
    coeffs[202] = d[195];
    coeffs[203] = d[201];
    coeffs[204] = d[210];
    coeffs[205] = d[197];
    coeffs[206] = d[203];
    coeffs[207] = d[212];
    coeffs[208] = d[205];
    coeffs[209] = d[214];
    coeffs[210] = d[198];
    coeffs[211] = d[204];
    coeffs[212] = d[213];
    coeffs[213] = d[206];
    coeffs[214] = d[215];
    coeffs[215] = d[217];
    coeffs[216] = d[207];
    coeffs[217] = d[216];
    coeffs[218] = d[218];
    coeffs[219] = d[219];
    coeffs[220] = d[220];
    coeffs[221] = d[224];
    coeffs[222] = d[232];
    coeffs[223] = d[244];
    coeffs[224] = d[221];
    coeffs[225] = d[225];
    coeffs[226] = d[233];
    coeffs[227] = d[245];
    coeffs[228] = d[228];
    coeffs[229] = d[236];
    coeffs[230] = d[248];
    coeffs[231] = d[222];
    coeffs[232] = d[226];
    coeffs[233] = d[234];
    coeffs[234] = d[246];
    coeffs[235] = d[229];
    coeffs[236] = d[237];
    coeffs[237] = d[249];
    coeffs[238] = d[240];
    coeffs[239] = d[252];
    coeffs[240] = d[223];
    coeffs[241] = d[227];
    coeffs[242] = d[235];
    coeffs[243] = d[247];
    coeffs[244] = d[230];
    coeffs[245] = d[238];
    coeffs[246] = d[250];
    coeffs[247] = d[241];
    coeffs[248] = d[253];
    coeffs[249] = d[256];
    coeffs[250] = d[231];
    coeffs[251] = d[239];
    coeffs[252] = d[251];
    coeffs[253] = d[242];
    coeffs[254] = d[254];
    coeffs[255] = d[257];
    coeffs[256] = d[243];
    coeffs[257] = d[255];
    coeffs[258] = d[258];
    coeffs[259] = d[259];
    coeffs[260] = d[260];
    coeffs[261] = d[264];
    coeffs[262] = d[272];
    coeffs[263] = d[284];
    coeffs[264] = d[261];
    coeffs[265] = d[265];
    coeffs[266] = d[273];
    coeffs[267] = d[285];
    coeffs[268] = d[268];
    coeffs[269] = d[276];
    coeffs[270] = d[288];
    coeffs[271] = d[262];
    coeffs[272] = d[266];
    coeffs[273] = d[274];
    coeffs[274] = d[286];
    coeffs[275] = d[269];
    coeffs[276] = d[277];
    coeffs[277] = d[289];
    coeffs[278] = d[280];
    coeffs[279] = d[292];
    coeffs[280] = d[263];
    coeffs[281] = d[267];
    coeffs[282] = d[275];
    coeffs[283] = d[287];
    coeffs[284] = d[270];
    coeffs[285] = d[278];
    coeffs[286] = d[290];
    coeffs[287] = d[281];
    coeffs[288] = d[293];
    coeffs[289] = d[296];
    coeffs[290] = d[271];
    coeffs[291] = d[279];
    coeffs[292] = d[291];
    coeffs[293] = d[282];
    coeffs[294] = d[294];
    coeffs[295] = d[297];
    coeffs[296] = d[283];
    coeffs[297] = d[295];
    coeffs[298] = d[298];
    coeffs[299] = d[299];
    coeffs[300] = d[300];
    coeffs[301] = d[303];
    coeffs[302] = d[309];
    coeffs[303] = d[318];
    coeffs[304] = d[301];
    coeffs[305] = d[304];
    coeffs[306] = d[310];
    coeffs[307] = d[319];
    coeffs[308] = d[306];
    coeffs[309] = d[312];
    coeffs[310] = d[321];
    coeffs[311] = d[302];
    coeffs[312] = d[305];
    coeffs[313] = d[311];
    coeffs[314] = d[320];
    coeffs[315] = d[307];
    coeffs[316] = d[313];
    coeffs[317] = d[322];
    coeffs[318] = d[315];
    coeffs[319] = d[324];
    coeffs[320] = d[308];
    coeffs[321] = d[314];
    coeffs[322] = d[323];
    coeffs[323] = d[316];
    coeffs[324] = d[325];
    coeffs[325] = d[327];
    coeffs[326] = d[317];
    coeffs[327] = d[326];
    coeffs[328] = d[328];
    coeffs[329] = d[329];
    coeffs[330] = d[330];
    coeffs[331] = d[333];
    coeffs[332] = d[339];
    coeffs[333] = d[348];
    coeffs[334] = d[331];
    coeffs[335] = d[334];
    coeffs[336] = d[340];
    coeffs[337] = d[349];
    coeffs[338] = d[336];
    coeffs[339] = d[342];
    coeffs[340] = d[351];
    coeffs[341] = d[332];
    coeffs[342] = d[335];
    coeffs[343] = d[341];
    coeffs[344] = d[350];
    coeffs[345] = d[337];
    coeffs[346] = d[343];
    coeffs[347] = d[352];
    coeffs[348] = d[345];
    coeffs[349] = d[354];
    coeffs[350] = d[338];
    coeffs[351] = d[344];
    coeffs[352] = d[353];
    coeffs[353] = d[346];
    coeffs[354] = d[355];
    coeffs[355] = d[357];
    coeffs[356] = d[347];
    coeffs[357] = d[356];
    coeffs[358] = d[358];
    coeffs[359] = d[359];



	// Setup elimination template
	static const int coeffs0_ind[] = { 0,40,80,110,150,190,220,260,300,330,1,41,81,111,151,191,221,261,301,331,4,44,80,84,114,154,190,194,224,264,300,304,330,334,10,50,90,120,160,200,230,270,310,340,5,45,81,85,115,155,191,195,225,265,301,305,331,335,11,51,84,91,121,161,194,201,231,271,304,311,334,341,17,57,90,97,127,167,200,207,237,277,310,317,340,347,2,42,82,112,152,192,222,262,302,332,3,43,83,113,153,193,223,263,303,333,8,48,88,118,158,198,228,268,308,338,9,49,89,119,159,199,229,269,309,339,18,58,98,128,168,208,238,278,318,348,19,59,99,129,169,209,239,279,319,349,29,69,105,139,179,215,249,289,325,355 };
	static const int coeffs1_ind[] = { 6,46,82,86,116,156,192,196,226,266,302,306,332,336,7,47,83,87,117,157,193,197,227,267,303,307,333,337,12,52,85,92,122,162,195,202,232,272,305,312,335,342,13,53,86,93,123,163,196,203,233,273,306,313,336,343,14,54,87,94,124,164,197,204,234,274,307,314,337,344,15,55,88,95,125,165,198,205,235,275,308,315,338,345,16,56,89,96,126,166,199,206,236,276,309,316,339,346,20,60,91,130,170,201,240,280,311,341,21,61,92,131,171,202,241,281,312,342,22,62,93,132,172,203,242,282,313,343,23,63,94,133,173,204,243,283,314,344,24,64,95,100,134,174,205,210,244,284,315,320,345,350,25,65,96,101,135,175,206,211,245,285,316,321,346,351,26,66,97,102,136,176,207,212,246,286,317,322,347,352,27,67,98,103,137,177,208,213,247,287,318,323,348,353,28,68,99,104,138,178,209,214,248,288,319,324,349,354,30,70,100,140,180,210,250,290,320,350,31,71,101,141,181,211,251,291,321,351,32,72,102,142,182,212,252,292,322,352,33,73,103,106,143,183,213,216,253,293,323,326,353,356,34,74,104,107,144,184,214,217,254,294,324,327,354,357,35,75,105,108,145,185,215,218,255,295,325,328,355,358,36,76,106,146,186,216,256,296,326,356,37,77,107,147,187,217,257,297,327,357,38,78,108,109,148,188,218,219,258,298,328,329,358,359,39,79,109,149,189,219,259,299,329,359 };
	static const int C0_ind[] = { 0,1,3,4,5,7,8,9,11,13,14,15,17,18,19,21,22,23,25,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,45,46,47,49,50,51,53,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,101,102,103,105,106,107,109,111,112,113,115,116,117,119,120,121,123,125,126,127,129,130,131,133,134,135,137,139,140,141,143,144,145,147,148,149,151,153,154,155,157,158,159,161,162,163,165,167,168,169,171,172,173,175,176,177,179,181,182,183,185,186,187,189,190,191,193,195 } ;
	static const int C1_ind[] = { 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,102,103,104,106,107,108,110,112,113,114,116,117,118,120,121,122,124,126,127,128,130,131,132,134,135,136,138,140,141,142,144,145,146,148,149,150,152,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,203,204,205,206,207,208,209,210,211,212,213,214,215,216,217,218,219,220,221,222,223,224,225,226,228,229,230,232,233,234,236,238,239,240,242,243,244,246,247,248,250,252,253,254,256,257,258,260,261,262,264,266,267,268,269,270,271,272,273,274,275,276,277,278,279,280,281,282,283,284,285,286,287,288,289,290,291,292,293,294,295,296,297,298,299,300,301,302,303,304,305,306,307,308,309,310,312,313,314,316,317,318,320,322,323,324,326,327,328,330,331,332,334,336,337,338,339,340,341,342,343,344,345,346,347,348,349,350,351,352,354,355,356,358,359,360,362 };

	Matrix<double,14,14> C0; C0.setZero();
	Matrix<double,14,26> C1; C1.setZero();
	for (int i = 0; i < 156; i++) { C0(C0_ind[i]) = coeffs(coeffs0_ind[i]); }
	for (int i = 0; i < 324; i++) { C1(C1_ind[i]) = coeffs(coeffs1_ind[i]); } 

	Matrix<double,14,26> C12 = C0.partialPivLu().solve(C1);




	// Setup action matrix
	Matrix<double,36, 26> RR;
	RR << -C12.bottomRows(10), Matrix<double,26,26>::Identity(26, 26);

	static const int AM_ind[] = { 3,4,0,10,11,5,6,1,12,13,14,15,16,2,7,8,21,22,23,24,25,9,29,30,31,34 };
	Matrix<double, 26, 26> AM;
	for (int i = 0; i < 26; i++) {
		AM.row(i) = RR.row(AM_ind[i]);
	}

	Matrix<std::complex<double>, 3, 26> sols;
	sols.setZero();

	// Solve eigenvalue problem


	double p[1+26];
	Matrix<double, 26, 26> AMp = AM;
	charpoly_danilevsky_piv(AMp, p);	
	double roots[26];
	int nroots;
	// find_real_roots_sturm(p, 26, roots, &nroots, 8, 0);
    nroots = realRoots (p, roots);
	fast_eigenvector_solver(roots, nroots, AM, sols);






	return sols;
}
// Action =  z
// Quotient ring basis (V) = x*y^2*z^2,y^3*z^2,x^2*y*z,x*y^2*z,y^3*z,x^2*z^2,x*y*z^2,x^3,x^2*y,x*y^2,y^3,x^2*z,x*y*z,y^2*z,x*z^2,y*z^2,x^2,x*y,y^2,x*z,y*z,z^2,x,y,z,1,
// Available monomials (RR*V) = x^2*y*z^2,x^3*z,y^2*z^2,x*y^2*z^3,y^3*z^3,x^2*z^3,x*y*z^3,x*z^3,y*z^3,z^3,x*y^2*z^2,y^3*z^2,x^2*y*z,x*y^2*z,y^3*z,x^2*z^2,x*y*z^2,x^3,x^2*y,x*y^2,y^3,x^2*z,x*y*z,y^2*z,x*z^2,y*z^2,x^2,x*y,y^2,x*z,y*z,z^2,x,y,z,1,


void fast_eigenvector_solver(double * eigv, int neig, Eigen::Matrix<double,26,26> &AM, Matrix<std::complex<double>,3,26> &sols) {
	static const int ind[] = { 0,1,2,5,6,7,13,14,15,21 };	
	// Truncated action matrix containing non-trivial rows
	Matrix<double, 10, 26> AMs;
	double zi[3];
	
	for (int i = 0; i < 10; i++)	{
		AMs.row(i) = AM.row(ind[i]);
	}
	for (int i = 0; i < neig; i++) {
		zi[0] = eigv[i];
		for (int j = 1; j < 3; j++)
		{
			zi[j] = zi[j - 1] * eigv[i];
		}
		Matrix<double, 10,10> AA;
        AA.col(0) = AMs.col(7);
        AA.col(1) = zi[0] * AMs.col(2) + AMs.col(8);
        AA.col(2) = zi[1] * AMs.col(5) + zi[0] * AMs.col(11) + AMs.col(16);
        AA.col(3) = zi[1] * AMs.col(0) + zi[0] * AMs.col(3) + AMs.col(9);
        AA.col(4) = zi[1] * AMs.col(6) + zi[0] * AMs.col(12) + AMs.col(17);
        AA.col(5) = zi[1] * AMs.col(14) + zi[0] * AMs.col(19) + AMs.col(22);
        AA.col(6) = zi[1] * AMs.col(1) + zi[0] * AMs.col(4) + AMs.col(10);
        AA.col(7) = zi[0] * AMs.col(13) + AMs.col(18);
        AA.col(8) = zi[1] * AMs.col(15) + zi[0] * AMs.col(20) + AMs.col(23);
        AA.col(9) = zi[1] * AMs.col(21) + zi[0] * AMs.col(24) + AMs.col(25);
        AA(5,0) = AA(5,0) - zi[0];
        AA(2,1) = AA(2,1) - zi[1];
        AA(3,2) = AA(3,2) - zi[2];
        AA(0,3) = AA(0,3) - zi[2];
        AA(4,4) = AA(4,4) - zi[2];
        AA(7,5) = AA(7,5) - zi[2];
        AA(1,6) = AA(1,6) - zi[2];
        AA(6,7) = AA(6,7) - zi[1];
        AA(8,8) = AA(8,8) - zi[2];
        AA(9,9) = AA(9,9) - zi[2];


		Matrix<double, 9, 1>  s = AA.leftCols(9).colPivHouseholderQr().solve(-AA.col(9));
        sols(1,i) = s(5);
        sols(2,i) = s(8);
        sols(0,i) = zi[0];

	}
}



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if (nrhs != 1) {
		mexErrMsgIdAndTxt("automatic_generator_cvpr:new_kE:nrhs", "One input required.");
	}
	if (nlhs != 1) {
		mexErrMsgIdAndTxt("automatic_generator_cvpr:new_kE:nlhs", "One output required.");
	}    
	if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0])) {
		mexErrMsgIdAndTxt("automatic_generator_cvpr:new_kE:notDouble", "Input data must be type double.");
	}
	if(mxGetNumberOfElements(prhs[0]) % 360 != 0) {
		mexErrMsgIdAndTxt("automatic_generator_cvpr:new_kE:incorrectSize", "Input size must be multiple of 360.");
	}
	int n_instances = mxGetNumberOfElements(prhs[0]) / 360;
	double *input = mxGetPr(prhs[0]);
	plhs[0] = mxCreateDoubleMatrix(3,26*n_instances,mxCOMPLEX);
	double* zr = mxGetPr(plhs[0]);
	double* zi = mxGetPi(plhs[0]);
	for(int k = 0; k < n_instances; k++) {
		const VectorXd data = Map<const VectorXd>(input + k*360, 360);
		MatrixXcd sols = solver_new_kE(data);
		Index offset = k*sols.size();
		for (Index i = 0; i < sols.size(); i++) {
			zr[i+offset] = sols(i).real();
			zi[i+offset] = sols(i).imag();
		}
	}
}


