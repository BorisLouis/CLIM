function F=exponential(x,data)

F = x(2) .* exp(-(data/x(1)))+x(3);