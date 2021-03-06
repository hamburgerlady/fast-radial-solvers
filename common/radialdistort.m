function y = radialdistort(x, kappa)
% y = radialdistort(x, kappa) applies radial distortion to the
% (homogeneous or inhomogeneous) coordinates x using the parameter kappa
% (should be negative) with the model y = (1 + kappa*|y|^2)*x. for
% coordinates in [-1,1]^2 kappa=-.01 is a mild distortion and kappa=-.5
% is a pretty heavy distortion.

ishom = (size(x, 1) == 3);
if(ishom)
    x = pflat(x);
    x = x(1:2, :);
end

% compute undistorted radius
ru2 = sum(x.^2);
ru = sqrt(ru2); 

% compute distorted radius
if(kappa == 0)
    rd = ru;
else
    beta = 1/4/kappa^2./ru2 - 1/kappa;
    beta = max(beta, 0); % to avoid sqrt() of a negative number in the next line.
    rd = 1/2/kappa./ru - sign(kappa)*sqrt(beta);
end

% compute distorted coordinates
y = repmat(rd./ru, size(x, 1), 1) .* x;

if(ishom)
    y = [y; ones(1, size(y, 2))];
end