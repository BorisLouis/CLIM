function C = contourc4(x,y,Z,v)
% Limited functionality wrapper to calculate output similar to contourc
% using bwperimtrace. Only acceptable input format is (x,y,Z,v) where
% x, y give coordinates for centre of pixels whose data are in 2D matrix Z,
% and v contains one or more contour values.
%
% Daniel Warren 03/02/2014
% Department of Oncology
% University of Oxford

v = unique(v);

for i = 1:length(v)
    bwpout{i} = bwperimtrace(Z >= v(i),[x(1) x(end)],[y(1) y(end)]);
end

C = [];
for i = 1:length(v)
    for j = 1:length(bwpout{i})
        C(end+1,1) = v(i);
        C(end,2) = size(bwpout{i}{j},1);
        C(end+(1:C(end,2)),:) = bwpout{i}{j};
    end
end

C = C';
 
end