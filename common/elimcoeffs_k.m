function [Ar,Br,Cr,Dr] =  elimcoeffs_k(A,B,C,D,nelim,ordo)


np = size(A,1);

A = A(:,ordo);
B = B(:,ordo);
C = C(:,ordo);
D = D(:,ordo);


AA = A(:,1:nelim);
AAs = A(1:nelim,1:nelim);
A = A(:,nelim+1:end);
B = B(:,nelim+1:end);
C = C(:,nelim+1:end);
D = D(:,nelim+1:end);

Ar = (AAs\A(1:nelim,:));
Br = (AAs\B(1:nelim,:));
Cr = (AAs\C(1:nelim,:));
Dr = (AAs\D(1:nelim,:));

for iii = nelim+1:np
    Ar = [Ar;A(iii,:)-AA(iii,1:nelim)*Ar(1:nelim,:)];
    Br = [Br;B(iii,:)-AA(iii,1:nelim)*Br(1:nelim,:)];
    Cr = [Cr;C(iii,:)-AA(iii,1:nelim)*Cr(1:nelim,:)];
    Dr = [Dr;D(iii,:)-AA(iii,1:nelim)*Dr(1:nelim,:)];
end

