
% where is your eigen?
cg_eigen_dir = '/usr/local/include/eigen3/'; 

% where is your sturm/charpoly ? 
template_dir = 'include/'; 

% output
outdir = 'solvers/';


%% kfEk
% which solver to compile?
fname = 'solvers/solver_new_kfEk.cpp';

mex(['-I"' cg_eigen_dir '"'],['-I"' template_dir '"'],'-outdir',outdir,'-O',fname)

%mex lincoeffs_kfE.cpp
%mex eqscoeffs_kfE.cpp
%mex Fsol_kfE.cpp

%% kE
% which solver to compile?
fname = 'solvers/solver_new_kE.cpp';

mex(['-I"' cg_eigen_dir '"'],['-I"' template_dir '"'],'-outdir',outdir,'-O',fname)

%% kfE
% which solver to compile?
fname = 'solvers/solver_new_kfE.cpp';

mex(['-I"' cg_eigen_dir '"'],['-I"' template_dir '"'],'-outdir',outdir,'-O',fname)

%% kEf
% which solver to compile?
fname = 'solvers/solver_new_kEf.cpp';

mex(['-I"' cg_eigen_dir '"'],['-I"' template_dir '"'],'-outdir',outdir,'-O',fname)

%% kfEf
% which solver to compile?
fname = 'solvers/solver_new_kfEf.cpp';

mex(['-I"' cg_eigen_dir '"'],['-I"' template_dir '"'],'-outdir',outdir,'-O',fname)


%% kfEfk
% which solver to compile?
fname = 'solvers/solver_new_kfEfk.cpp';

mex(['-I"' cg_eigen_dir '"'],['-I"' template_dir '"'],'-outdir',outdir,'-O',fname)


%% kF
% which solver to compile?
fname = 'solvers/solver_new_kF.cpp';

mex(['-I"' cg_eigen_dir '"'],['-I"' template_dir '"'],'-outdir',outdir,'-O',fname)

%% kFk
% which solver to compile?
fname = 'solvers/solver_new_kFk.cpp';

mex(['-I"' cg_eigen_dir '"'],['-I"' template_dir '"'],'-outdir',outdir,'-O',fname)



%% k2fEk1
% which solver to compile?
fname = 'solvers/solver_new_k2fEk1.cpp';

mex(['-I"' cg_eigen_dir '"'],['-I"' template_dir '"'],'-outdir',outdir,'-O',fname)

%% k2fEfk1
% which solver to compile?
fname = 'solvers/solver_new_k2fEfk1.cpp';

mex(['-I"' cg_eigen_dir '"'],['-I"' template_dir '"'],'-outdir',outdir,'-O',fname)

%% k2Fk1
% which solver to compile?
fname = 'solvers/solver_new_k2Fk1.cpp';

mex(['-I"' cg_eigen_dir '"'],['-I"' template_dir '"'],'-outdir',outdir,'-O',fname)



