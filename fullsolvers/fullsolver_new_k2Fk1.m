function [F,k2,k1] = fullsolver_new_k2Fk1(u,v)
% function [F,k2,k1] = fullsolver_new_k2Fk1(u,v)
%
% input: image coordinates 2x9 u (first image) and 2x9 v (second image)
% output: Estimated fundamental matrices 9xn F and radial 1xn k2 and k1 (n solutions)
% 
% Magnus Oskarsson 2021

ktype = 3;
nelim = 4;
ordo = [1 2 4 5 7 8 3 6 9];

sz = [ 81 16];
[A,B,C,D] = lincoeffs_k(u,v,ktype);
[Ar,Br,Cr,Dr] = elimcoeffs_k(A,B,C,D,nelim,ordo);

Br = Br(:,3:end);
Cr = Cr(:,[1 2 5]);
Dr = Dr(:,end);
data = data_from_lin_k2Fk1_cc(Ar,Br,Cr,Dr);
data = normalize_data_eqs(data,sz);
sols = solver_new_k2Fk1(data);
sols(:,sum(abs(sols))<1e-10)=[];

k2 = sols(1,:);
k1 = sols(2,:);

F = F_from_sol_k2Fk1_cc(Ar(:),Br(:),Cr(:),Dr(:),k2(:),k1(:));

