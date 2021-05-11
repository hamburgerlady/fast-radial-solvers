%% init
nrsam = 10000;
binny = -16:.3:3;
testcases = cell(0);

%% testcase kE
caseinfo.name = 'kE';
caseinfo.nvars = 3;
caseinfo.np = 6;
caseinfo.ktype = 1;
caseinfo.ftype = 0;
caseinfo.fullsolver = str2func('fullsolver_new_kE');

testcases{end+1} = caseinfo;


%% testcase kfE
caseinfo.name = 'kfE';
caseinfo.sz = [  24    30    35    35];
caseinfo.nvars = 2;
caseinfo.np = 7;
caseinfo.ordo = [1 2 4 5 7 8 3 6 9];
caseinfo.ktype = 1;
caseinfo.ftype = 1;
caseinfo.nelim = 6;
caseinfo.datafun = str2func('data_from_lin_kfE_cc');
caseinfo.solver = str2func('solver_new_kfE');

testcases{end+1} = caseinfo;

%% testcase kEf
caseinfo.name = 'kEf';
caseinfo.sz = [  24    40    40    45];
caseinfo.nvars = 2;
caseinfo.np = 7;
caseinfo.ordo = [1 2 4 5 7 8 3 6 9];
caseinfo.ktype = 1;
caseinfo.ftype = 4;
caseinfo.nelim = 6;
caseinfo.datafun = str2func('data_from_lin_kEf_cc');
caseinfo.solver = str2func('solver_new_kEf');

testcases{end+1} = caseinfo;

%% testcase kfEf
caseinfo.name = 'kfEf';
caseinfo.sz = [24 60];
caseinfo.nvars = 2;
caseinfo.np = 7;
caseinfo.ordo = [1 2 4 5 7 8 3 6 9];
caseinfo.ktype = 1;
caseinfo.ftype = 2;
caseinfo.nelim = 6;
caseinfo.datafun = str2func('data_from_lin_kfEf_cc');
caseinfo.solver = str2func('solver_new_kfEf');

testcases{end+1} = caseinfo;

%% testcase kfEk
caseinfo.name = 'kfEk';
caseinfo.sz = [50 65 80 85];
caseinfo.nvars = 2;
caseinfo.np = 7;
caseinfo.ordo = [1 2 4 5 7 8 3 6 9];
caseinfo.ktype = 2;
caseinfo.ftype = 1;
caseinfo.nelim = 4;
caseinfo.datafun = str2func('data_from_lin_kfEk_cc');
caseinfo.solver = str2func('solver_new_kfEk');

testcases{end+1} = caseinfo;



%% testcase kfEfk
caseinfo.name = 'kfEfk';
caseinfo.sz = [50 117];
caseinfo.nvars = 2;
caseinfo.np = 7;
caseinfo.ordo = [1 2 4 5 7 8 3 6 9];
caseinfo.ktype = 2;
caseinfo.ftype = 2;
caseinfo.nelim = 4;
caseinfo.datafun = str2func('data_from_lin_kfEfk_cc');
caseinfo.solver = str2func('solver_new_kfEfk');

testcases{end+1} = caseinfo;

%% testcase k2Ek1
caseinfo.name = 'k2Ek1';
caseinfo.sz = [ 238   238   204   238   238   204   210   210   180   180     7];
caseinfo.nvars = 3;
caseinfo.np = 7;
caseinfo.ordo = [1 2 4 5 7 8 3 6 9];
caseinfo.ktype = 3;
caseinfo.ftype = 0;
caseinfo.nelim = 4;
caseinfo.datafun = str2func('data_from_lin_k2Ek1_cc');
caseinfo.solver = str2func('solver_new_k2Ek1');
testcases{end+1} = caseinfo;


%% testcase kF
caseinfo.name = 'kF';
caseinfo.sz = 9;
caseinfo.nvars = 1;
caseinfo.np = 8;
caseinfo.ordo = [1 2 4 5 7 8 3 6 9];
caseinfo.ktype = 1;
caseinfo.ftype = 3;
caseinfo.nelim = 6;
caseinfo.datafun = str2func('data_from_lin_kF_cc');
caseinfo.solver = str2func('solver_new_kF');

testcases{end+1} = caseinfo;



%% testcase kFk
caseinfo.name = 'kFk';
caseinfo.sz = 17;
caseinfo.nvars = 1;
caseinfo.np = 8;
caseinfo.ordo = [1 2 4 5 7 8 3 6 9];
caseinfo.ktype = 2;
caseinfo.ftype = 3;
caseinfo.nelim = 4;
caseinfo.datafun = str2func('data_from_lin_kFk_cc');
caseinfo.solver = str2func('solver_new_kFk');

testcases{end+1} = caseinfo;



%% testcase k2fEk1
caseinfo.name = 'k2fEk1';
caseinfo.sz = [ 81   132   132   143];
caseinfo.nvars = 2;
caseinfo.np = 8;
caseinfo.ordo = [1 2 4 5 7 8 3 6 9];
caseinfo.ktype = 3;
caseinfo.ftype = 1;
caseinfo.nelim = 4;
caseinfo.datafun = str2func('data_from_lin_k2fEk1_cc');
caseinfo.solver = str2func('solver_new_k2fEk1');

testcases{end+1} = caseinfo;


%% testcase k2fEfk1
caseinfo.name = 'k2fEfk1';
caseinfo.sz = [ 81   221];
caseinfo.nvars = 2;
caseinfo.np = 8;
caseinfo.ordo = [1 2 4 5 7 8 3 6 9];
caseinfo.ktype = 3;
caseinfo.ftype = 2;
caseinfo.nelim = 4;
caseinfo.datafun = str2func('data_from_lin_k2fEfk1_cc');
caseinfo.solver = str2func('solver_new_k2fEfk1');

testcases{end+1} = caseinfo;

%% testcase k2Fk1
caseinfo.name = 'k2Fk1';
caseinfo.sz = [81 16];
caseinfo.nvars = 2;
caseinfo.np = 9;
caseinfo.ordo = [1 2 4 5 7 8 3 6 9];
caseinfo.ktype = 3;
caseinfo.ftype = 3;
caseinfo.nelim = 4;
caseinfo.datafun = str2func('data_from_lin_k2Fk1_cc');
caseinfo.solver = str2func('solver_new_k2Fk1');
testcases{end+1} = caseinfo;


%% test cases

nrc = length(testcases);
allnormF = nan(nrc,nrsam);

for cc = 1:nrc
    disp([cc nrc])
    for iii = 1:nrsam
        
        [u,v,gt,K1,K2,P1,P2,E,F] = getRandomData(caseinfo.ktype,caseinfo.ftype,caseinfo.nvars,caseinfo.np);
        Fsol = caseinfo.fullsolver(u,v);
        
        if ~isempty(Fsol)
            allnormF(cc,iii)  = min(sqrt(sum((Fsol-F(:)).^2)));
        end
    end
    
    figure(cc)
    clf
        histogram(log10(abs(allnormF(cc,:))),binny,'Normalization','pdf')
        
    title(caseinfo.name)
    
end






