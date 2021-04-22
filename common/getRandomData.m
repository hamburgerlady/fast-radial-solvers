function [u,v,gt,K1,K2,P1,P2,E,F] = getRandomData(ktype,ftype,nvars,np)


switch ftype
    case 0 % E
        f2 = 1;
        f1 = 1;
    case 1 % fE       
        f2 = 1+rand;
        f1 = 1;
    case 2 % fEf
        f2 = 1+rand;
        f1 = f2;
    case 3 % f2Ef1 = F
        f2 = 1+rand;
        f1 = 1+rand;
    case 4 % Ef
        f2 = 1;
        f1 = 1+rand;
end

switch ktype
    case 1 % kF
        k2 = -rand/10;
        k1 = 0;
    case 2 % kFk
        k2 =  -rand/10;
        k1 = k2;
    case 3 % k2Fk1
        k2 = -rand/10;
        k1 = -rand/10;
end


K1 = diag([f1 f1 1]);
K2 = diag([f2 f2 1]);

[R,~]=qr(randn(3));
t = randn(3,1)*10;
P1 = K1*[eye(3) zeros(3,1)];
P2 = K2*[R t];

E = crossm(t)*R;
F = K2*E*K1;
F = F/F(end);

x1 = zeros(3,np);
x2 = zeros(3,np);
found = 0;
while found<np
    X = [randn(3,1)*10;1];
    x1(:,found+1)=P1*X;
    x2(:,found+1)=P2*X;
    if x1(3,found+1)>0 && x2(3,found+1)>0
        found = found+1;
    end
end
x1 = x1./x1(3,:);
x2 = x2./x2(3,:);
y1 = radialdistort(x1,k1);
y2 = radialdistort(x2,k2);

u = y1(1:2,:);
v = y2(1:2,:);

if ktype == 3
    gt = [k2 k1];
else
    gt = k2;
end
fsmall = F([6 3]);
gt = [gt fsmall(1:(nvars-length(gt)))];
gt = gt(:);




