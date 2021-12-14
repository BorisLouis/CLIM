function contourcplot(C,col)
% Plot contourc output.
%
% Daniel Warren 03/02/2014
% Department of Oncology
% University of Oxford

if ~exist('col','var')
    col = 'b';
end

i = 1;
hold on;
while i < size(C,2)
    plot(C(1,i+(1:C(2,i))),C(2,i+(1:C(2,i))),col);
    i = i + C(2,i) + 1;
end

clabel(C);
hold off;

end