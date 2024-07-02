function [FitPar,Fit]=expFit(A,domain)

A=A(:);
domain =domain(:);
[val,ind] = max(A);

tauGuess = prctile(domain,20);
AGuess = max(A);
bkgGuess = 0;


expEq = fittype('a*exp(-(x/t))+b');
startPoints = [AGuess bkgGuess tauGuess];
options = fitoptions(expEq);
options.StartPoint = startPoints;
options.Lower = [0.5,0, domain(2),];
options.Upper = [1,0.5,domain(end)];
options.Exclude = A<0;

weight = ones(length(A),1);
weight(domain>20) = 0.5;
options.weight = weight;

f1 = fit(domain,A,expEq,options);

Fit = f1(domain);
FitPar = [f1.t,f1.a, f1.b];

% %                   tau                       A             y0       y0           
% lb        = [abs(domain(1)-domain(2))     min(A)           min(A)];
% ub        = [domain(end)     max(A)+0.1*max(abs(A))      val];
% initguess = [      tauGuess                 AGuess                         bkgGuess];
% opts = optimset('Display','off');
% FitPar=lsqcurvefit(@SimpleFitting.exponential,initguess,domain,A,lb,ub,opts);
% 
% Fit= SimpleFitting.exponential(FitPar,domain);
