function [F,k] = fullsolver_new_kfEf(u,v)
% function [F,k] = fullsolver_new_kfEf(u,v)
%
% input: image coordinates 2x7 u (first image) and 2x7 v (second image)
% output: Estimated fundamental matrices 9xn F and radial 1xn k (n solutions)


ktype = 1;
nelim = 6;
ordo = [1 2 4 5 7 8 3 6 9];

sz = [ 24 60];
[A,B,C,D] = lincoeffs_k(u,v,ktype);
[Ar,Br,~,~] = elimcoeffs_k(A,B,C,D,nelim,ordo);
data = data_from_lin_kfEf_cc(Ar,Br);
data = normalize_data_eqs(data,sz);
sols = solver_new_kfEf(data);
sols(:,sum(abs(sols))<1e-10)=[];

k = sols(1,:);
f1 = sols(2,:);
F = F_from_sol_kfEf_cc(Ar(:),Br(:),k(:),f1(:));

