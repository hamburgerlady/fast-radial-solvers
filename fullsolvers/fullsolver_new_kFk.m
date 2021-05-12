function [F,k] = fullsolver_new_kFk(u,v)
% function [F,k] = fullsolver_new_kFk(u,v)
%
% input: image coordinates 2x8 u (first image) and 2x8 v (second image)
% output: Estimated fundamental matrices 9xn F and radial 1xn k (n solutions)
% 
% Magnus Oskarsson 2021

ktype = 2;
nelim = 4;
ordo = [1 2 4 5 7 8 3 6 9];


sz = 17;
[A,B,C,D] = lincoeffs_k(u,v,ktype);
[Ar,Br,~,Dr] = elimcoeffs_k(A,B,C,D,nelim,ordo);
Dr = Dr(:,end);
data = data_from_lin_kFk_cc(Ar,Br,Dr);
data = normalize_data_eqs(data,sz);
sols = solver_new_kFk(data);
sols((abs(sols))<1e-10)=[];

k = sols(1,:);

F = F_from_sol_kFk_cc(Ar(:),Br(:),Dr(:),k(:));

