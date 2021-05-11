function [F,k] = fullsolver_new_kF(u,v)
% function [F,k] = fullsolver_new_kF(u,v)
%
% input: image coordinates 2x8 u (first image) and 2x8 v (second image)
% output: Estimated fundamental matrices 9xn F and radial 1xn k (n solutions)

ktype = 1;
nelim = 6;
ordo = [1 2 4 5 7 8 3 6 9];

sz = 9;
[A,B,C,D] = lincoeffs_k(u,v,ktype);
[Ar,Br,~,~] = elimcoeffs_k(A,B,C,D,nelim,ordo);
data = data_from_lin_kF_cc(Ar,Br);
data = normalize_data_eqs(data,sz);
sols = solver_new_kF(data);
sols((abs(sols))<1e-10)=[];

k = sols(1,:);

F = F_from_sol_kF_cc(Ar(:),Br(:),k(:));

