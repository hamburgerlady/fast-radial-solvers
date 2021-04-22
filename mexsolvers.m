
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
