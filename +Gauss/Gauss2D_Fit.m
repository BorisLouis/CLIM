%Fitting
function [gpar, fit]=Gauss2D_Fit(A,domain)

radius= 1.5*(size(A,2));
Bkg0 = 0;
wguess = 3;

BkgUB=0.9*max(A(:));

[~,ind] = max(A(:));
[idy0,idx0] = ind2sub(size(A),ind);
x0 = domain(1,idx0,1);
y0 = domain(idy0,1,2);

domX = domain(:,:,1);
domY = domain(:,:,2);

minX = min(domX(:));
minY = min(domY(:));

maxX = max(domX(:));
maxY = max(domY(:));

bkglb=0;

PSF=1; %pixel

curvefitoptions = optimset('Display','off'); 

            %[peak Int,     sigma X,    sigma Y,        Bkg,
            
lb        = [0              PSF           PSF               bkglb ... 
    minX    minY    0 ]; 
%mean X,         mean Y angle

ub        = [1.5*max(A(:))  radius      radius          BkgUB ...
    maxX    maxY    pi];
initguess = [max(A(:))      wguess      wguess          Bkg0  ...
    x0              y0              0];


gpar=lsqcurvefit(@Gauss2DRot,initguess,domain,A,lb,ub,curvefitoptions);

fit = Gauss2DRot(gpar,domain);