function [FitPar,Fit]=gauss1D(A,domain,guess)

switch nargin
    case 2 
        sigGuess = abs((abs(domain(2))-abs(domain(1))))*3;

        [val,idx] = max(A);
        muGuess = domain(idx);
    case 3
        sigGuess = guess.sig;
        muGuess  = guess.mu;
    otherwise
        error('wrong number of arguments');
end

A=A(:);
domain =domain(:);
[val,ind] = max(A);

%                   Sigma                       mu             A        y0           
lb        = [0     0     0        0];
ub        = [abs(domain(1)-domain(2))*20     max(domain)+0.1*max(abs(domain))     3*val      val];
initguess = [      sigGuess                  muGuess                          val-min(A) min(A)];
opts = optimset('Display','off');
FitPar=lsqcurvefit(@gaussian,initguess,domain,A,lb,ub,opts);

Fit=FitPar(3)*exp(-((domain-FitPar(2))./(sqrt(2).*FitPar(1))).^2)+FitPar(4);
