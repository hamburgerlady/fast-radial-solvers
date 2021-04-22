function [A,B,C,D] = lincoeffs_k(u,v,ktype)

% u,v: 2 x np
% A,B,C,D : np x 9
% A:1 B:k1 C:k2 D:k1k2
if nargin<3
    ktype = 1;
end
np = size(u,2);
%u2 = sum(u.^2);
%v2 = sum(v.^2);
A = zeros(np,9);
B = zeros(np,9);
C = zeros(np,9);
D = zeros(np,9);

for iii = 1:np
    u1 = u(1,iii);
    u2 = u(2,iii);
    v1 = v(1,iii);
    v2 = v(2,iii);
    switch ktype
        case 1 % kF
            A(iii,:) = [u1*v1, u1*v2, u1, u2*v1, u2*v2, u2, v1, v2, 1];
            B(iii,:) = [0, 0, u1*(v1^2 + v2^2), 0, 0, u2*(v1^2 + v2^2), 0, 0, v1^2 + v2^2];
        case 2 % kFk
            A(iii,:) = [u1*v1, u1*v2, u1, u2*v1, u2*v2, u2, v1, v2, 1];
            B(iii,:) = [0, 0, u1*(v1^2 + v2^2), 0, 0, u2*(v1^2 + v2^2), v1*(u1^2 + u2^2), v2*(u1^2 + u2^2), u1^2 + u2^2 + v1^2 + v2^2];
            D(iii,:) = [0, 0, 0, 0, 0, 0, 0, 0, (u1^2 + u2^2)*(v1^2 + v2^2)];
        case 3 % k1Fk2
            A(iii,:) = [u1*v1, u1*v2, u1, u2*v1, u2*v2, u2, v1, v2, 1];
            B(iii,:) = [0, 0, u1*(v1^2 + v2^2), 0, 0, u2*(v1^2 + v2^2), 0, 0, v1^2 + v2^2];
            C(iii,:) = [0, 0, 0, 0, 0, 0, v1*(u1^2 + u2^2), v2*(u1^2 + u2^2), u1^2 + u2^2];
            D(iii,:) = [0, 0, 0, 0, 0, 0, 0, 0, (u1^2 + u2^2)*(v1^2 + v2^2)];
    end
end

            
            
    
    