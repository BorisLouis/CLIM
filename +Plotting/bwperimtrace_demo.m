function bwperimtrace_demo
% Demonstration of the bwperimtrace function. Creates and displays a random
% posterized image, feeds each level into bwperimtrace and then shows the
% output boundaries both separately and as an overlay.
%
% Daniel Warren 05/12/2015
% Department of Oncology
% University of Oxford

nlevels = 12;
nfeatures = 10;
xsize = 100;
ysize = 250;

xv = (1:xsize)-(xsize/2);
yv = (1:ysize)-(ysize/2);
[x y] = meshgrid(xv,yv);

im = zeros(size(x));
for i = 1:nfeatures
    im = im+(-1+2*mod(i,2))*rand*sqrt((x+xsize*(rand-0.5)).^2+(y+ysize*(rand-0.5)).^2);
end
im = round((nlevels-1)*(im-min(im(:)))/(max(im(:))-min(im(:))));

bounds = cell(0);
for t = unique(im(:))'
    bounds{end+1} = bwperimtrace(im == t,[xv(1) xv(end)],[yv(1) yv(end)]);
end

colors = jet(nlevels);

figure;
subplot(1,3,1);
imshow(im,[0 nlevels-1],'Xdata',[xv(1) xv(end)],'Ydata',[yv(1) yv(end)]);
subplot(1,3,2);
hold on;
for i = 1:nlevels
    for j = 1:length(bounds{i})
        plot(bounds{i}{j}(:,1),bounds{i}{j}(:,2),'Color',colors(i,:))
    end
end
axis image ij;
hold off;
subplot(1,3,3);
hold on;
imshow(im,[0 nlevels-1],'Xdata',[xv(1) xv(end)],'Ydata',[yv(1) yv(end)]);
for i = 1:nlevels
    for j = 1:length(bounds{i})
        plot(bounds{i}{j}(:,1),bounds{i}{j}(:,2),'Color',colors(i,:))
    end
end
axis image;
hold off;

end