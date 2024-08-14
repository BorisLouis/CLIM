
%% Compare samples

% Look at size distribution
% Intensity distribution
% Sil distribution?
pxSize = 0.2;
pxArea = pxSize^2;



%% Size distribution O2
%bigGrain

figure
nameO2 = {'BigO2','PorO2','MedO2','SmallO2'};
nameN2 = {'BigN2','PorN2','MedN2','SmallN2'};

for i = 1: length(nameO2)

    lb = min([(nameO2{i}).results.nPx]*pxArea);
    ub = max([BigO2.results.nPx]*pxArea);
    edgesO2 = lb:0.2:ub;
    lb = min([BigN2.results.nPx]*pxArea);
    ub = max([BigN2.results.nPx]*pxArea);

    edgesN2 = lb:0.2:ub;

    [NO2,edgesO2] = histcounts([BigO2.results.nPx]*pxArea,edgesO2);
    [NN2,edgesN2] = histcounts([BigN2.results.nPx]*pxArea,edgesN2);

    subplot(1,4,1)
    plot(edgesO2(1:end-1),NO2./sum(NO2))
    hold on
    plot(edgesN2(1:end-1),NN2./sum(NN2))

end

