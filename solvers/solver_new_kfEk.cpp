#include <Eigen/Dense>
#include "mex.h"

#define DEG  66
#include "sturm_mart.h"
#include "charpoly.h"



using namespace Eigen;

void fast_eigenvector_solver(double * eigv, int neig, Eigen::MatrixXd &AM, Matrix<std::complex<double>,2,66> &sols);



MatrixXcd solver_new_kfEk(const VectorXd& data)
{
	// Compute coefficients
    const double* d = data.data();
    VectorXd coeffs(280);
    coeffs[0] = d[0];
    coeffs[1] = d[11];
    coeffs[2] = d[23];
    coeffs[3] = d[36];
    coeffs[4] = d[1];
    coeffs[5] = d[12];
    coeffs[6] = d[24];
    coeffs[7] = d[37];
    coeffs[8] = d[2];
    coeffs[9] = d[13];
    coeffs[10] = d[25];
    coeffs[11] = d[38];
    coeffs[12] = d[3];
    coeffs[13] = d[14];
    coeffs[14] = d[26];
    coeffs[15] = d[39];
    coeffs[16] = d[4];
    coeffs[17] = d[15];
    coeffs[18] = d[27];
    coeffs[19] = d[40];
    coeffs[20] = d[5];
    coeffs[21] = d[16];
    coeffs[22] = d[28];
    coeffs[23] = d[41];
    coeffs[24] = d[6];
    coeffs[25] = d[17];
    coeffs[26] = d[29];
    coeffs[27] = d[42];
    coeffs[28] = d[7];
    coeffs[29] = d[18];
    coeffs[30] = d[30];
    coeffs[31] = d[43];
    coeffs[32] = d[8];
    coeffs[33] = d[19];
    coeffs[34] = d[31];
    coeffs[35] = d[44];
    coeffs[36] = d[9];
    coeffs[37] = d[20];
    coeffs[38] = d[32];
    coeffs[39] = d[45];
    coeffs[40] = d[10];
    coeffs[41] = d[21];
    coeffs[42] = d[33];
    coeffs[43] = d[46];
    coeffs[44] = d[22];
    coeffs[45] = d[34];
    coeffs[46] = d[47];
    coeffs[47] = d[35];
    coeffs[48] = d[48];
    coeffs[49] = d[49];
    coeffs[50] = d[50];
    coeffs[51] = d[61];
    coeffs[52] = d[73];
    coeffs[53] = d[86];
    coeffs[54] = d[100];
    coeffs[55] = d[51];
    coeffs[56] = d[62];
    coeffs[57] = d[74];
    coeffs[58] = d[87];
    coeffs[59] = d[101];
    coeffs[60] = d[52];
    coeffs[61] = d[63];
    coeffs[62] = d[75];
    coeffs[63] = d[88];
    coeffs[64] = d[102];
    coeffs[65] = d[53];
    coeffs[66] = d[64];
    coeffs[67] = d[76];
    coeffs[68] = d[89];
    coeffs[69] = d[103];
    coeffs[70] = d[54];
    coeffs[71] = d[65];
    coeffs[72] = d[77];
    coeffs[73] = d[90];
    coeffs[74] = d[104];
    coeffs[75] = d[55];
    coeffs[76] = d[66];
    coeffs[77] = d[78];
    coeffs[78] = d[91];
    coeffs[79] = d[105];
    coeffs[80] = d[56];
    coeffs[81] = d[67];
    coeffs[82] = d[79];
    coeffs[83] = d[92];
    coeffs[84] = d[106];
    coeffs[85] = d[57];
    coeffs[86] = d[68];
    coeffs[87] = d[80];
    coeffs[88] = d[93];
    coeffs[89] = d[107];
    coeffs[90] = d[58];
    coeffs[91] = d[69];
    coeffs[92] = d[81];
    coeffs[93] = d[94];
    coeffs[94] = d[108];
    coeffs[95] = d[59];
    coeffs[96] = d[70];
    coeffs[97] = d[82];
    coeffs[98] = d[95];
    coeffs[99] = d[109];
    coeffs[100] = d[60];
    coeffs[101] = d[71];
    coeffs[102] = d[83];
    coeffs[103] = d[96];
    coeffs[104] = d[110];
    coeffs[105] = d[72];
    coeffs[106] = d[84];
    coeffs[107] = d[97];
    coeffs[108] = d[111];
    coeffs[109] = d[85];
    coeffs[110] = d[98];
    coeffs[111] = d[112];
    coeffs[112] = d[99];
    coeffs[113] = d[113];
    coeffs[114] = d[114];
    coeffs[115] = d[115];
    coeffs[116] = d[129];
    coeffs[117] = d[144];
    coeffs[118] = d[160];
    coeffs[119] = d[177];
    coeffs[120] = d[116];
    coeffs[121] = d[130];
    coeffs[122] = d[145];
    coeffs[123] = d[161];
    coeffs[124] = d[178];
    coeffs[125] = d[117];
    coeffs[126] = d[131];
    coeffs[127] = d[146];
    coeffs[128] = d[162];
    coeffs[129] = d[179];
    coeffs[130] = d[118];
    coeffs[131] = d[132];
    coeffs[132] = d[147];
    coeffs[133] = d[163];
    coeffs[134] = d[180];
    coeffs[135] = d[119];
    coeffs[136] = d[133];
    coeffs[137] = d[148];
    coeffs[138] = d[164];
    coeffs[139] = d[181];
    coeffs[140] = d[120];
    coeffs[141] = d[134];
    coeffs[142] = d[149];
    coeffs[143] = d[165];
    coeffs[144] = d[182];
    coeffs[145] = d[121];
    coeffs[146] = d[135];
    coeffs[147] = d[150];
    coeffs[148] = d[166];
    coeffs[149] = d[183];
    coeffs[150] = d[122];
    coeffs[151] = d[136];
    coeffs[152] = d[151];
    coeffs[153] = d[167];
    coeffs[154] = d[184];
    coeffs[155] = d[123];
    coeffs[156] = d[137];
    coeffs[157] = d[152];
    coeffs[158] = d[168];
    coeffs[159] = d[185];
    coeffs[160] = d[124];
    coeffs[161] = d[138];
    coeffs[162] = d[153];
    coeffs[163] = d[169];
    coeffs[164] = d[186];
    coeffs[165] = d[125];
    coeffs[166] = d[139];
    coeffs[167] = d[154];
    coeffs[168] = d[170];
    coeffs[169] = d[187];
    coeffs[170] = d[126];
    coeffs[171] = d[140];
    coeffs[172] = d[155];
    coeffs[173] = d[171];
    coeffs[174] = d[188];
    coeffs[175] = d[127];
    coeffs[176] = d[141];
    coeffs[177] = d[156];
    coeffs[178] = d[172];
    coeffs[179] = d[189];
    coeffs[180] = d[128];
    coeffs[181] = d[142];
    coeffs[182] = d[157];
    coeffs[183] = d[173];
    coeffs[184] = d[190];
    coeffs[185] = d[143];
    coeffs[186] = d[158];
    coeffs[187] = d[174];
    coeffs[188] = d[191];
    coeffs[189] = d[159];
    coeffs[190] = d[175];
    coeffs[191] = d[192];
    coeffs[192] = d[176];
    coeffs[193] = d[193];
    coeffs[194] = d[194];
    coeffs[195] = d[195];
    coeffs[196] = d[210];
    coeffs[197] = d[226];
    coeffs[198] = d[243];
    coeffs[199] = d[261];
    coeffs[200] = d[196];
    coeffs[201] = d[211];
    coeffs[202] = d[227];
    coeffs[203] = d[244];
    coeffs[204] = d[262];
    coeffs[205] = d[197];
    coeffs[206] = d[212];
    coeffs[207] = d[228];
    coeffs[208] = d[245];
    coeffs[209] = d[263];
    coeffs[210] = d[198];
    coeffs[211] = d[213];
    coeffs[212] = d[229];
    coeffs[213] = d[246];
    coeffs[214] = d[264];
    coeffs[215] = d[199];
    coeffs[216] = d[214];
    coeffs[217] = d[230];
    coeffs[218] = d[247];
    coeffs[219] = d[265];
    coeffs[220] = d[200];
    coeffs[221] = d[215];
    coeffs[222] = d[231];
    coeffs[223] = d[248];
    coeffs[224] = d[266];
    coeffs[225] = d[201];
    coeffs[226] = d[216];
    coeffs[227] = d[232];
    coeffs[228] = d[249];
    coeffs[229] = d[267];
    coeffs[230] = d[202];
    coeffs[231] = d[217];
    coeffs[232] = d[233];
    coeffs[233] = d[250];
    coeffs[234] = d[268];
    coeffs[235] = d[203];
    coeffs[236] = d[218];
    coeffs[237] = d[234];
    coeffs[238] = d[251];
    coeffs[239] = d[269];
    coeffs[240] = d[204];
    coeffs[241] = d[219];
    coeffs[242] = d[235];
    coeffs[243] = d[252];
    coeffs[244] = d[270];
    coeffs[245] = d[205];
    coeffs[246] = d[220];
    coeffs[247] = d[236];
    coeffs[248] = d[253];
    coeffs[249] = d[271];
    coeffs[250] = d[206];
    coeffs[251] = d[221];
    coeffs[252] = d[237];
    coeffs[253] = d[254];
    coeffs[254] = d[272];
    coeffs[255] = d[207];
    coeffs[256] = d[222];
    coeffs[257] = d[238];
    coeffs[258] = d[255];
    coeffs[259] = d[273];
    coeffs[260] = d[208];
    coeffs[261] = d[223];
    coeffs[262] = d[239];
    coeffs[263] = d[256];
    coeffs[264] = d[274];
    coeffs[265] = d[209];
    coeffs[266] = d[224];
    coeffs[267] = d[240];
    coeffs[268] = d[257];
    coeffs[269] = d[275];
    coeffs[270] = d[225];
    coeffs[271] = d[241];
    coeffs[272] = d[258];
    coeffs[273] = d[276];
    coeffs[274] = d[242];
    coeffs[275] = d[259];
    coeffs[276] = d[277];
    coeffs[277] = d[260];
    coeffs[278] = d[278];
    coeffs[279] = d[279];



	// Setup elimination template
	static const int coeffs0_ind[] = { 0,1,51,116,196,1,2,52,117,197,2,3,53,118,198,3,54,119,199,5,6,57,122,1,2,52,202,117,6,7,58,123,2,3,53,203,118,7,59,124,3,54,204,119,10,11,63,128,6,7,58,3,53,208,123,2,11,64,129,7,59,54,209,124,3,15,69,134,11,64,59,54,214,129,7,3,19,74,139,15,69,64,59,3,54,219,134,11,7,23,79,144,19,74,69,64,7,59,3,224,139,15,11,180,100,95,40,260,175,36,100,265,180,40,0,50,115,195,4,5,56,121,0,1,51,201,116,9,10,62,127,5,6,57,2,52,207,122,1,14,15,68,133,10,11,63,7,58,53,3,213,128,6,2,27,84,149,23,79,74,69,11,64,7,229,144,19,15 };
	static const int coeffs1_ind[] = { 4,55,120,0,50,200,115,8,60,125,4,55,0,50,205,120,8,9,61,126,4,5,56,1,51,206,121,0,12,65,130,8,60,4,55,50,0,210,125,12,13,66,131,8,9,61,5,56,51,1,211,126,4,0,13,14,67,132,9,10,62,6,57,52,2,212,127,5,1,16,70,135,12,65,8,60,55,50,4,215,130,0,16,17,71,136,12,13,66,9,61,56,0,51,5,216,131,8,4,1,17,18,72,137,13,14,67,10,62,57,1,52,6,217,132,9,5,2,18,19,73,138,14,15,68,11,63,58,2,53,7,218,133,10,6,3,20,75,140,16,70,12,65,60,55,8,220,135,4,20,21,76,141,16,17,71,13,66,61,4,56,0,9,221,136,12,8,5,21,22,77,142,17,18,72,14,67,62,5,57,1,10,222,137,13,9,6,22,23,78,143,18,19,73,15,68,63,6,58,2,11,223,138,14,10,7,24,80,145,20,75,16,70,65,60,12,225,140,8,24,25,81,146,20,21,76,17,71,66,8,61,4,13,226,141,16,12,9,25,26,82,147,21,22,77,18,72,67,9,62,5,14,227,142,17,13,10,26,27,83,148,22,23,78,19,73,68,10,63,6,15,228,143,18,14,11,28,85,150,24,80,20,75,70,65,16,230,145,12,28,29,86,151,24,25,81,21,76,71,12,66,8,17,231,146,20,16,13,29,30,87,152,25,26,82,22,77,72,13,67,9,18,232,147,21,17,14,30,31,88,153,26,27,83,23,78,73,14,68,10,19,233,148,22,18,15,31,89,154,27,84,79,74,15,69,11,234,149,23,19,32,90,155,28,85,24,80,75,70,20,235,150,16,32,33,91,156,28,29,86,25,81,76,16,71,12,21,236,151,24,20,17,33,34,92,157,29,30,87,26,82,77,17,72,13,22,237,152,25,21,18,34,35,93,158,30,31,88,27,83,78,18,73,14,23,238,153,26,22,19,35,94,159,31,89,84,79,19,74,15,239,154,27,23,36,95,160,32,90,28,85,80,75,24,240,155,20,36,37,96,161,32,33,91,29,86,81,20,76,16,25,241,156,28,24,21,37,38,97,162,33,34,92,30,87,82,21,77,17,26,242,157,29,25,22,38,39,98,163,34,35,93,31,88,83,22,78,18,27,243,158,30,26,23,39,99,164,35,94,89,84,23,79,19,244,159,31,27,40,100,165,36,95,32,90,85,80,28,245,160,24,40,41,101,166,36,37,96,33,91,86,24,81,20,29,246,161,32,28,25,41,42,102,167,37,38,97,34,92,87,25,82,21,30,247,162,33,29,26,42,43,103,168,38,39,98,35,93,88,26,83,22,31,248,163,34,30,27,43,104,169,39,99,94,89,27,84,23,249,164,35,31,170,40,100,36,95,90,85,32,250,165,28,44,105,171,40,41,101,37,96,91,28,86,24,33,251,166,36,32,29,44,45,106,172,41,42,102,38,97,92,29,87,25,34,252,167,37,33,30,45,46,107,173,42,43,103,39,98,93,30,88,26,35,253,168,38,34,31,46,108,174,43,104,99,94,31,89,27,254,169,39,35,175,40,100,95,90,36,255,170,32,176,44,105,41,101,96,32,91,28,37,256,171,40,36,33,47,109,177,44,45,106,42,102,97,33,92,29,38,257,172,41,37,34,47,48,110,178,45,46,107,43,103,98,34,93,30,39,258,173,42,38,35,48,111,179,46,108,104,99,35,94,31,259,174,43,39,181,44,105,101,36,96,32,41,261,176,40,37,182,47,109,45,106,102,37,97,33,42,262,177,44,41,38,49,112,183,47,48,110,46,107,103,38,98,34,43,263,178,45,42,39,49,113,184,48,111,108,104,39,99,35,264,179,46,43,185,105,40,101,36,44,266,181,41,186,47,109,106,41,102,37,45,267,182,44,42,187,49,112,48,110,107,42,103,38,46,268,183,47,45,43,114,188,49,113,111,108,43,104,39,269,184,48,46,105,40,270,185,44,189,109,44,106,41,47,271,186,45,190,49,112,110,45,107,42,48,272,187,47,46,191,114,113,111,46,108,43,273,188,49,48,109,44,274,189,47,192,112,47,110,45,49,275,190,48,193,114,113,48,111,46,276,191,49,112,47,277,192,49,194,114,49,113,48,278,193,114,49,279,194 };
	static const int C0_ind[] = { 0,1,2,3,14,19,20,21,22,33,38,39,40,41,52,57,59,60,71,76,77,78,79,80,81,82,90,91,95,96,97,98,99,100,101,109,110,114,116,117,118,120,128,129,133,134,135,136,137,138,139,140,141,147,148,149,152,154,155,156,158,160,166,167,168,171,173,174,175,177,179,180,185,186,187,188,190,192,193,194,196,198,199,200,201,204,205,206,207,209,211,212,213,215,217,218,219,220,221,223,224,225,226,231,237,239,241,242,243,246,258,261,262,265,267,268,269,280,285,286,287,288,289,290,291,299,300,304,305,306,307,308,309,310,311,312,318,319,320,323,324,325,326,327,328,329,330,331,332,336,337,338,339,340,342,344,345,346,348,350,351,352,353,354,356,357,358,359 } ;
	static const int C1_ind[] = { 1,2,3,5,6,14,15,20,21,22,24,25,26,27,33,34,38,39,40,41,42,43,44,45,46,52,53,54,58,59,60,62,63,64,65,66,70,71,72,76,77,78,79,80,81,82,83,84,85,89,90,91,92,93,95,96,97,98,99,100,101,102,103,104,108,109,110,111,112,115,116,117,119,120,121,122,123,125,127,128,129,132,133,134,135,136,137,138,139,140,141,142,143,144,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,184,185,186,187,188,189,191,192,193,195,196,197,198,199,201,203,204,205,208,209,210,211,212,213,214,215,216,217,218,219,220,221,222,223,224,225,226,227,228,229,230,231,232,233,234,235,236,237,238,239,240,241,242,243,244,245,246,247,248,249,250,251,252,253,254,255,256,257,258,259,260,261,262,263,264,265,267,268,269,271,272,273,274,275,277,279,280,281,284,285,286,287,288,289,290,291,292,293,294,295,296,297,298,299,300,301,302,303,304,305,306,307,308,309,310,311,312,313,314,315,316,317,318,319,320,321,322,323,324,325,326,327,328,329,330,331,332,333,334,335,336,337,338,339,340,341,343,344,345,347,348,349,350,351,353,355,356,357,360,361,362,363,364,365,366,367,368,369,370,371,372,373,374,375,376,377,378,379,380,381,382,383,384,385,386,387,388,389,390,391,392,393,394,395,396,397,398,399,400,401,402,403,404,405,406,407,408,409,410,411,412,413,414,415,416,417,418,420,421,422,424,426,427,428,429,430,432,433,434,435,438,439,440,442,443,444,445,446,448,450,451,452,455,456,457,458,459,460,461,462,463,464,465,466,467,468,469,470,471,472,473,474,475,476,477,478,479,480,481,482,483,484,485,486,487,488,489,490,491,492,493,494,495,496,497,498,499,500,501,502,503,504,505,506,507,508,509,510,511,512,513,515,516,517,519,521,522,523,524,525,527,528,529,530,533,534,535,537,538,539,540,541,543,545,546,547,550,551,552,553,554,555,556,557,558,559,560,561,562,563,564,565,566,567,568,569,570,571,572,573,574,575,576,577,578,579,580,581,582,583,584,585,586,587,588,589,590,591,592,593,594,595,596,597,598,599,600,601,602,603,604,605,606,607,608,610,611,612,614,616,617,618,619,620,622,623,624,625,628,629,630,632,633,634,635,636,638,640,641,642,645,646,647,648,649,650,651,652,653,654,655,656,657,658,659,660,661,662,663,664,665,666,667,668,669,670,671,672,673,674,675,676,677,678,679,680,681,682,683,684,685,686,687,688,689,690,691,692,693,694,695,696,697,698,699,700,701,702,703,705,706,707,709,711,712,713,714,715,717,718,719,720,725,727,728,729,730,731,733,735,736,737,740,742,743,744,745,746,747,748,749,750,751,752,753,754,755,756,757,758,759,760,761,762,763,764,765,766,767,768,769,770,771,772,773,774,775,776,777,778,779,780,781,782,783,784,785,786,787,788,789,790,791,792,793,794,795,796,797,798,800,801,802,804,806,807,808,809,810,812,813,814,815,820,824,825,826,828,830,831,832,835,839,841,842,843,844,845,846,847,848,849,850,851,852,853,854,856,857,858,859,860,861,862,863,864,865,866,867,868,869,870,871,872,873,874,875,876,877,878,879,880,881,882,883,884,885,886,887,888,889,890,891,892,893,895,896,897,899,901,902,903,904,905,907,908,909,910,915,919,920,921,922,923,924,925,926,927,929,930,934,936,937,938,939,940,941,942,943,944,945,946,947,948,949,951,952,953,954,955,956,957,958,959,960,961,962,963,964,965,966,967,968,969,971,972,973,975,977,978,979,980,981,983,984,985,986,991,997,998,999,1000,1001,1002,1003,1006,1010,1014,1015,1016,1017,1018,1019,1020,1021,1022,1024,1025,1029,1031,1032,1033,1034,1035,1036,1037,1038,1039,1040,1041,1042,1043,1044,1047,1048,1049,1051,1053,1054,1055,1056,1057,1059,1060,1061,1062,1075,1076,1078,1079,1082,1086,1092,1093,1094,1095,1096,1097,1098,1101,1105,1109,1110,1111,1112,1113,1114,1115,1116,1117,1119,1120,1124,1127,1129,1130,1131,1132,1133,1135,1136,1137,1138,1151,1152,1154,1155,1158,1162,1168,1169,1170,1171,1172,1173,1174,1177,1181,1186,1187,1188,1189,1190,1192,1193,1195,1208,1209,1211,1212,1215,1219,1225,1226,1227,1228,1230,1231,1246,1247,1249,1250 };

	Matrix<double,19,19> C0; C0.setZero();
	Matrix<double,19,66> C1; C1.setZero();
	for (int i = 0; i < 168; i++) { C0(C0_ind[i]) = coeffs(coeffs0_ind[i]); }
	for (int i = 0; i < 952; i++) { C1(C1_ind[i]) = coeffs(coeffs1_ind[i]); } 

	Matrix<double,19,66> C12 = C0.partialPivLu().solve(C1);




	// Setup action matrix
	// Matrix<double,71, 66> RR;
    MatrixXd RR(71, 66);	
	RR << -C12.bottomRows(5), Matrix<double,66,66>::Identity(66, 66);

	static const int AM_ind[] = { 0,5,1,6,7,2,8,9,10,3,11,12,13,14,15,16,17,18,19,20,21,22,4,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,49,50,51,52,53,54,55,56,57,58,59,60,62,63,64,66,67,69 };
	// Matrix<double, 66, 66> AM;
    MatrixXd AM(66, 66);
	for (int i = 0; i < 66; i++) {
		AM.row(i) = RR.row(AM_ind[i]);
	}

	Matrix<std::complex<double>, 2, 66> sols;
	sols.setZero();

	// Solve eigenvalue problem


	double p[1+66];
	// Matrix<double, 66, 66> AMp = AM;
    MatrixXd AMp = AM;
	charpoly_danilevsky_piv(AMp, p);	
	double roots[66];
	int nroots;
	// find_real_roots_sturm(p, 66, roots, &nroots, 8, 0);
    nroots = realRoots (p, roots);
	fast_eigenvector_solver(roots, nroots, AM, sols);






	return sols;
}
// Action =  y
// Quotient ring basis (V) = x^4*y^13,x^4*y^12,x^3*y^13,x^4*y^11,x^3*y^12,x^2*y^13,x^4*y^10,x^3*y^11,x^2*y^12,x*y^13,x^4*y^9,x^3*y^10,x^2*y^11,x*y^12,x^4*y^8,x^3*y^9,x^2*y^10,x*y^11,x^4*y^7,x^3*y^8,x^2*y^9,x*y^10,y^11,x^4*y^6,x^3*y^7,x^2*y^8,x*y^9,y^10,x^4*y^5,x^3*y^6,x^2*y^7,x*y^8,y^9,x^4*y^4,x^3*y^5,x^2*y^6,x*y^7,y^8,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^4*y^2,x^3*y^3,x^2*y^4,x*y^5,y^6,x^3*y^2,x^2*y^3,x*y^4,y^5,x^3*y,x^2*y^2,x*y^3,y^4,x^3,x^2*y,x*y^2,y^3,x^2,x*y,y^2,x,y,1,
// Available monomials (RR*V) = x^4*y^14,x^3*y^14,x^2*y^14,x*y^14,y^12,x^4*y^13,x^4*y^12,x^3*y^13,x^4*y^11,x^3*y^12,x^2*y^13,x^4*y^10,x^3*y^11,x^2*y^12,x*y^13,x^4*y^9,x^3*y^10,x^2*y^11,x*y^12,x^4*y^8,x^3*y^9,x^2*y^10,x*y^11,x^4*y^7,x^3*y^8,x^2*y^9,x*y^10,y^11,x^4*y^6,x^3*y^7,x^2*y^8,x*y^9,y^10,x^4*y^5,x^3*y^6,x^2*y^7,x*y^8,y^9,x^4*y^4,x^3*y^5,x^2*y^6,x*y^7,y^8,x^4*y^3,x^3*y^4,x^2*y^5,x*y^6,y^7,x^4*y^2,x^3*y^3,x^2*y^4,x*y^5,y^6,x^3*y^2,x^2*y^3,x*y^4,y^5,x^3*y,x^2*y^2,x*y^3,y^4,x^3,x^2*y,x*y^2,y^3,x^2,x*y,y^2,x,y,1,


void fast_eigenvector_solver(double * eigv, int neig, Eigen::MatrixXd &AM, Matrix<std::complex<double>,2,66> &sols) {
	static const int ind[] = { 0,2,5,9,22 };	
	// Truncated action matrix containing non-trivial rows
	Matrix<double, 5, 66> AMs;
	double zi[14];
	
	for (int i = 0; i < 5; i++)	{
		AMs.row(i) = AM.row(ind[i]);
	}
	for (int i = 0; i < neig; i++) {
		zi[0] = eigv[i];
		for (int j = 1; j < 14; j++)
		{
			zi[j] = zi[j - 1] * eigv[i];
		}
		Matrix<double, 5,5> AA;
        AA.col(0) = zi[12] * AMs.col(0) + zi[11] * AMs.col(1) + zi[10] * AMs.col(3) + zi[9] * AMs.col(6) + zi[8] * AMs.col(10) + zi[7] * AMs.col(14) + zi[6] * AMs.col(18) + zi[5] * AMs.col(23) + zi[4] * AMs.col(28) + zi[3] * AMs.col(33) + zi[2] * AMs.col(38) + zi[1] * AMs.col(43);
        AA.col(1) = zi[12] * AMs.col(2) + zi[11] * AMs.col(4) + zi[10] * AMs.col(7) + zi[9] * AMs.col(11) + zi[8] * AMs.col(15) + zi[7] * AMs.col(19) + zi[6] * AMs.col(24) + zi[5] * AMs.col(29) + zi[4] * AMs.col(34) + zi[3] * AMs.col(39) + zi[2] * AMs.col(44) + zi[1] * AMs.col(48) + zi[0] * AMs.col(52) + AMs.col(56);
        AA.col(2) = zi[12] * AMs.col(5) + zi[11] * AMs.col(8) + zi[10] * AMs.col(12) + zi[9] * AMs.col(16) + zi[8] * AMs.col(20) + zi[7] * AMs.col(25) + zi[6] * AMs.col(30) + zi[5] * AMs.col(35) + zi[4] * AMs.col(40) + zi[3] * AMs.col(45) + zi[2] * AMs.col(49) + zi[1] * AMs.col(53) + zi[0] * AMs.col(57) + AMs.col(60);
        AA.col(3) = zi[12] * AMs.col(9) + zi[11] * AMs.col(13) + zi[10] * AMs.col(17) + zi[9] * AMs.col(21) + zi[8] * AMs.col(26) + zi[7] * AMs.col(31) + zi[6] * AMs.col(36) + zi[5] * AMs.col(41) + zi[4] * AMs.col(46) + zi[3] * AMs.col(50) + zi[2] * AMs.col(54) + zi[1] * AMs.col(58) + zi[0] * AMs.col(61) + AMs.col(63);
        AA.col(4) = zi[10] * AMs.col(22) + zi[9] * AMs.col(27) + zi[8] * AMs.col(32) + zi[7] * AMs.col(37) + zi[6] * AMs.col(42) + zi[5] * AMs.col(47) + zi[4] * AMs.col(51) + zi[3] * AMs.col(55) + zi[2] * AMs.col(59) + zi[1] * AMs.col(62) + zi[0] * AMs.col(64) + AMs.col(65);
        AA(0,0) = AA(0,0) - zi[13];
        AA(1,1) = AA(1,1) - zi[13];
        AA(2,2) = AA(2,2) - zi[13];
        AA(3,3) = AA(3,3) - zi[13];
        AA(4,4) = AA(4,4) - zi[11];


		Matrix<double, 4, 1>  s = AA.leftCols(4).colPivHouseholderQr().solve(-AA.col(4));
        sols(1,i) = s(3);
        sols(0,i) = zi[0];

	}
}



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if (nrhs != 1) {
		mexErrMsgIdAndTxt("automatic_generator_cvpr:new_kfEk:nrhs", "One input required.");
	}
	if (nlhs != 1) {
		mexErrMsgIdAndTxt("automatic_generator_cvpr:new_kfEk:nlhs", "One output required.");
	}    
	if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0])) {
		mexErrMsgIdAndTxt("automatic_generator_cvpr:new_kfEk:notDouble", "Input data must be type double.");
	}
	if(mxGetNumberOfElements(prhs[0]) % 280 != 0) {
		mexErrMsgIdAndTxt("automatic_generator_cvpr:new_kfEk:incorrectSize", "Input size must be multiple of 280.");
	}
	int n_instances = mxGetNumberOfElements(prhs[0]) / 280;
	double *input = mxGetPr(prhs[0]);
	plhs[0] = mxCreateDoubleMatrix(2,66*n_instances,mxCOMPLEX);
	double* zr = mxGetPr(plhs[0]);
	double* zi = mxGetPi(plhs[0]);
	for(int k = 0; k < n_instances; k++) {
		const VectorXd data = Map<const VectorXd>(input + k*280, 280);
		MatrixXcd sols = solver_new_kfEk(data);
		Index offset = k*sols.size();
		for (Index i = 0; i < sols.size(); i++) {
			zr[i+offset] = sols(i).real();
			zi[i+offset] = sols(i).imag();
		}
	}
}


