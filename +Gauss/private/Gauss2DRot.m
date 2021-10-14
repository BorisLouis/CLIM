%Generate a gaussian

function F = Gauss2DRot(x,data)

X = data(:,:,1);
Y = data(:,:,2);

a =  cos(x(7))^2/(2*x(2)^2) + sin(x(7))^2/(2*x(3)^2);
b = -sin(2*x(7))/(4*x(2)^2) + sin(2*x(7))/(4*x(3)^2);
c =  sin(x(7))^2/(2*x(2)^2) + cos(x(7))^2/(2*x(3)^2);

F = x(1)*exp( - (a*(X-x(5)).^2 - 2*b*(X-x(5)).*(Y-x(6)) +...
    c*(Y-x(6)).^2)) + x(4);