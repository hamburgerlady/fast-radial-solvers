function [E,k] = fullsolver_new_kE(u,v)
% function [E,k] = fullsolver_new_kE(u,v)
%
% input: image coordinates 2x6 u (first image) and 2x6 v (second image)
% output: Estimated essential matrices 9xn E and radial 1xn k (n solutions)
% 
% Magnus Oskarsson 2021


ktype = 1;
nelim = 6;
ordo = [1 2 4 5 7 8 3 6 9];
iordo = [1 2 7 3 4 8 5 6 9];
sz = [40    40    30    40    40    30    40    40    30    30];
[A,B,C,D] = lincoeffs_k(u,v,ktype);
[Ar,Br,~,~] = elimcoeffs_k(A,B,C,D,nelim,ordo);
data = data_from_lin_kE_cc(Ar,Br);
data = normalize_data_eqs(data,sz);
sols = solver_new_kE(data);
sols(:,sum(abs(sols))<1e-10)=[];

vv = [sols(3,:);sols(2,:);ones(1,size(sols,2))];
k = sols(1,:);

Ev = -Ar*vv-Br*(vv.*k);
E = [Ev;vv];
E = E(iordo,:);
