%% Simulate correlated data

%% Simulation input
sizeIm = 64;
nFrames = 100;

nParticles = 4;
corrThreshold = 0.7;%smaller is more selective here (0 is perfect correlation)
model.name = 'gaussian';
model.sigma_x = 3;
model.sigma_y = 3;
r = 2; %radius for checking correlation


%% Simulations
data = zeros(sizeIm,sizeIm,nFrames);
coord = [28,28;36,36;28,36;36,28];
%coord = [32,32;96,96;96,32;32,96];
[X,Y] = meshgrid(1:sizeIm,1:sizeIm);
%coord = zeros(nParticles,2);
for i = 1:nParticles
    
    x0 = coord(i,1);
    y0 = coord(i,2);
    c  = 0;
    BaseInt = 3e4;
    secondInt = 2*BaseInt;
    int = BaseInt;
    for j = 1 : 100
        
        num = rand(1);
        if num > 0.7
            if int==BaseInt
                int=secondInt;
            else
                int = BaseInt;
            end
        else
        end
        
        PSF = Sim.getPSF(X,Y,x0,y0,model);
        sPSF = Sim.samplePSF(PSF,int,false);
        
        data(:,:,j) = data(:,:,j) + double(sPSF); 


    end
    coord(i,:) = [x0,y0];   
end

noise = randn(size(data));

finalData = data + ones(size(data))*100 +noise*20;

%% Clustering
[corrRel,corrSum]  = corrAnalysis.getCorrRelation2(finalData,r,corrThreshold);


%%
listCorrPx = reshape(corrRel,size(corrRel,1)*size(corrRel,2),1);
inds    = (1:length(listCorrPx))';

idx2Delete = cellfun(@isempty,listCorrPx);

listCorrPx(idx2Delete) =[];

inds(idx2Delete) = [];

[n,p] = ind2sub(size(corrRel),inds);
pxIntList = zeros(length(n),size(finalData,3));
for i =1:length(n)
   
    pxIntList(i,:) = finalData(n(i),p(i),:);
    
end

% get distance map
dist = pdist2([n,p],[n,p]);
%normalize the distance scale
dist = dist./max(dist(:));
distanceMap = 1-corrcoef(pxIntList');

%#1 Try clustering the data together.
testMask = zeros(size(finalData,1),size(finalData,2));
T1 = clusterdata(distanceMap,4);
for i = 1: length(T1)
    testMask(inds(i)) = T1(i);
    
end

figure
imagesc(testMask);

% 
% %#2 DBScan is very similar to what I implemented as pseudo clustering so it
% %fails for similar reasons.
% idx = dbscan(distanceMap,24,50,'Distance','precomputed');
% unique(idx)

