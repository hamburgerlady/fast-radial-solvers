function [F,k] = fullsolver_new_kfEk(u,v)
% function [F,k] = fullsolver_new_kfEk(u,v)
%
% input: image coordinates 2x7 u (first image) and 2x7 v (second image)
% output: Estimated fundamental matrices 9xn F and radial 1xn k (n solutions)


ktype = 2;
nelim = 4;
ordo = [1 2 4 5 7 8 3 6 9];

sz = [50 65 80 85];
[A,B,C,D] = lincoeffs_k(u,v,ktype);
[Ar,Br,~,Dr] = elimcoeffs_k(A,B,C,D,nelim,ordo);
Dr = Dr(:,end);
data = data_from_lin_kfEk_cc(Ar,Br,Dr);
data = normalize_data_eqs(data,sz);
sols = solver_new_kfEk(data);
sols(:,sum(abs(sols))<1e-10)=[];

k = sols(1,:);
f1 = sols(2,:);
F = F_from_sol_kfEk_cc(Ar(:),Br(:),Dr(:),k(:),f1(:));

