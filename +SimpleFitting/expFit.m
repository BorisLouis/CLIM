function [FitPar,Fit]=expFit(A,domain)

A=A(:);
domain =domain(:);
[val,ind] = max(A);

tauGuess = prctile(domain,20);
AGuess = max(A)-min(A);
bkgGuess = min(A);

%                   tau                       A             y0       y0           
lb        = [abs(domain(1)-domain(2))     min(A)           min(A)];
ub        = [domain(end)     max(A)+0.1*max(abs(A))      val];
initguess = [      tauGuess                 AGuess                         bkgGuess];
opts = optimset('Display','off');
FitPar=lsqcurvefit(@SimpleFitting.exponential,initguess,domain,A,lb,ub,opts);

Fit= SimpleFitting.exponential(FitPar,domain);
