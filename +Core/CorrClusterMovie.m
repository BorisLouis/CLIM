classdef CorrClusterMovie < Core.Movie
    %General definition of a movie (master class from which many other
    %movie type will inherit
    %The Movie hold information and path to the movie but does not store
    %any data inside. Methods allow to display and do stuff with the data.
    
    properties (SetAccess = 'private')
        corrMask
        nCluster
        driftCorr
    end
        
    methods
        function obj = CorrClusterMovie(raw,info,driftCorr)
            
            obj = obj@Core.Movie(raw,info);
            
            if nargin>2
                obj.driftCorr = driftCorr;
            else
                obj.driftCorr = true;
            end
                        
        end
        
        function correctDrift(obj)
            if obj.driftCorr
                disp('Correcting drift')
                disp('Check if drift corrected data exist')
                %#1 Check if drift was already corrected
                [run] = obj.checkDrift;

                if run 
                    correlationInfo.corrSz = 75; %in px. Radius of the ROI used for correlation
                    %correlation function
                    correlationInfo.driftPeriod = 10; %in Frame, Number of frame that are averaged
                    correlationInfo.maxDrift = 20;
                    %for driftCalculation ==> 1 mean that drift is calculated for each frame
                    scalingFactor = 1;%Used for interpolation in sub-pixel Drift correction 
                    %objects of interest
                    
                    disp('Loading Data ...')
                    h = waitbar(0,'Loading Frames');
                    nFrame = obj.raw.movInfo.maxFrame;
                    allFrames = zeros(obj.raw.movInfo.Length,obj.raw.movInfo.Width,nFrame);
                    for i = 1:nFrame

                        allFrames(:,:,i) = uint16(obj.getFrame(i));
                        waitbar(i/(nFrame),h,...
                            ['Loading Frames ' num2str(i) '/' num2str(nFrame)]);

                    end
                    close(h);
                    disp('Correcting drift...')
                    %fix Drift
                    [corrData,~] = PreProcess.CorrelationDrift(allFrames,scalingFactor,correlationInfo);

                    currentFolder = obj.raw.movInfo.Path;

                    newDir = [currentFolder filesep 'DriftCorr'];
                    mkdir(newDir);

                    fileName = [newDir filesep 'corrData.mat'];
                    corrData = uint16(corrData);
                    disp('Saving corrected data')
                    %dataStorage.nBTiff(fileName,uint16(corrData),16);
                    save(fileName,'corrData','-v7.3');
                else
                    disp('Corrected data found ! ')
                end
            else
                
                
                
            end
        end
        
        function [data] = getFrame(obj,idx)
            switch nargin
                case 1
                    
                    idx = 1:obj.raw.maxFrame(1);
                
                case 2

                    [idx] = Core.Movie.checkFrame(idx,obj.raw.maxFrame(1));
                otherwise
                    error('too many input arguments');
            end
            %if we drift corrected and if the data exist
            if and(obj.driftCorr,~obj.checkDrift)
                
                path = [obj.raw.movInfo.Path filesep 'DriftCorr' filesep 'corrData.mat'];
                assert(isfile(path),'No drift corrected data found, please run correctDrift first');
                
                dataObj = matfile(path);
                [data] = dataObj.corrData(:,:,idx);
                
            else
                
                [data] = getFrame@Core.Movie(obj,idx);
                
            end
            
        end
        
          
        function [corrData] = loadFrames(obj,frames)
            %simple method to load all requested frames and allow to take
            %roi
            assert(length(frames)>=100,'Frames requested is lower than 100 please process at least 100 fr');
            fr = obj.checkFrame(frames,obj.raw.movInfo.maxFrame);
                
            h = questdlg('Do you want to use a ROI?','Question to user','Yes','No', 'Yes');
            frame = obj.getFrame(1);
            if strcmp(h,'Yes')
                figure
                
                imagesc(frame(:,:,1))
                test = drawrectangle();
                ROI  = round(test.Position);    
            else
                ROI = [];
                
            end
            
            frames = obj.getFrame(fr);
             %  h = waitbar(0,'Loading Frames');
             %  frames = zeros(obj.raw.movInfo.Length,obj.raw.movInfo.Width,length(fr));
%             for i = fr(1):fr(end)
%                 
%                 frames(:,:,i) = obj.getFrame(i);
%                 waitbar((i-fr(1))/(fr(end)-fr(1)),h,['Loading Frames ' num2str(i) '/' num2str(length(fr))]);
%             end
%             close(h);
                    
            if isempty(ROI)
               ROI = [1 1 size(frames,2),size(frames,1)];
            end
            
            %Fix ROI
            if mod(ROI(3),2) ==0
                ROI(3) = ROI(3)-1;
            end
            
            if mod(ROI(4),2) ==0
                ROI(4) = ROI(4)-1;
            end
            %apply ROI
            corrData = frames(ROI(2):ROI(2)+ROI(4),ROI(1):ROI(1) + ROI(3),:);
       
        end
        
        function [corrMask] = getCorrelationMask(obj,data,corrInfo)
            r = corrInfo.r;
            corrThreshold = corrInfo.thresh;
            h=waitbar(0,'Normalizing Data...');
            % #1 Normalize the data
            data2Cluster = data-min(data,[],3);
            data2Cluster = data2Cluster./max(data2Cluster,[],3);
              
            waitbar(0.1,h,'Getting correlation relationships');
            %#2 Get correlation relationship between pixels
            [corrRel]  = obj.getCorrRelation(data2Cluster,r,corrThreshold);
            
            waitbar(0.4,h,'Preparing data');
            %#3 reshap pixel relationship and get corresponding indices        
            listCorrPx = reshape(corrRel,size(corrRel,1)*size(corrRel,2),1);
            inds    = (1:length(listCorrPx))';
            
            %#4 Clean data by keeping only pixel that have correlation
            %relation
            idx2Delete = cellfun(@isempty,listCorrPx);
            listCorrPx(idx2Delete) =[];
            inds(idx2Delete) = [];
            
            %calculate distancemap
            [n,p] = ind2sub(size(corrRel),inds);
            pxIntList = zeros(length(n),size(data2Cluster,3));
            for i =1:length(n)

                pxIntList(i,:) = data2Cluster(n(i),p(i),:);

            end
            distanceMap = 1-corrcoef(pxIntList');
            waitbar(0.5,h,'Intensity fluctuation based clustering');
            %perform pseudo-clustering
            dim = size(data2Cluster);
            [corrMask] = obj.corrClustering(listCorrPx,inds,distanceMap,dim);

            figure
            imagesc(corrMask)
            axis image
            
            obj.corrMask = corrMask;
            obj.nCluster = max(corrMask(:));
            
        end
        
    end
    methods(Static)
        function [corrRel] = getCorrRelation(data2Cluster,r,corrThreshold)
            %function to find correlation relation between a each pixel of
            %an image and its neighbor pixels
            
            corrRel  = cell(size(data2Cluster,1),size(data2Cluster,2));
            %loop through pixels
            for i = 1:size(data2Cluster,1)
                for j = 1:size(data2Cluster,2)
                    %take a pixel
                    currentPxCoord = [i,j];
                    %find its neighbor
                    neighbor = Core.CorrClusterMovie.findNeighbor(currentPxCoord,size(data2Cluster),r);

                    corr = zeros(size(neighbor,1),1);                 
                    data1 = squeeze(data2Cluster(currentPxCoord(1),currentPxCoord(2),:));
                    %calculate correlation 1-pearson coefficient ==> 0 is
                    %correlated 2 is anti correated 1 is uncorrelated
                    for k = 1:size(neighbor,1)
                        data2 = squeeze(data2Cluster(neighbor(k,1),neighbor(k,2),:));
                        tmpCorr = corrcoef(data1,data2);
                        corr(k) = 1-tmpCorr(2,1);

                    end
                    %if correlation between pixel is sufficent we keep
                    %track of those pixel as being correlated to the
                    %current pixel
                    [neighborIdx] = sub2ind(size(corrRel),neighbor(:,1),neighbor(:,2));
                    currPxIdx     = sub2ind(size(corrRel),currentPxCoord(:,1),currentPxCoord(:,2));
                    corr(neighborIdx==currPxIdx) = [];
                    neighborIdx(neighborIdx==currPxIdx) = [];
                    
                    if all(corr>corrThreshold)

                    else
                        idx = neighborIdx(corr<corrThreshold,:);
                        
                        if ~isempty(idx)
                           %convert to indices for simplicity later
                           corrRel{currPxIdx} = idx;
                        end
                    end
                end
            end
        end
        
        
        function [corrMask] = corrClustering(listCorrPx,inds,distanceMap,dim)
            %function actually does the clustering.
            safeCount = 2*dim(1)*dim(2);
            %initialize the variable
            count = 1;
            group = 1;
            corrMask = zeros(dim(1),dim(2));
            treatedIdx = NaN(length(inds),1);
            indsCopy = inds;
            %as long as the list of pixel that are correlated is not empty
            %we keep going
            while ~isempty(listCorrPx) 
                %take the index of the first pixel to be treated
                currIndex = inds(1);
                %add to a new list the pixel that are correlated with the
                %currently treated pixel
                currList = listCorrPx{1};
                %check that the list does not contain already treated
                %pixels
                currList(ismember(currList,treatedIdx,'row'),:) =[];
                %if the list is empty then the pixel is "dead" and we
                %remove it from the list
                if isempty(currList)
                    listCorrPx(1) = [];
                    idx = find(isnan(treatedIdx(:,1)),1);
                    treatedIdx(idx) = currIndex;
                    inds(1) = [];
                else
                    %otherwise we treat all the pixel in the list in the
                    %same way
                    while ~isempty(currList)
                        % add the new pixel to the group and keep track of it
                        corrMask(currIndex) = group;
                        idx = find(isnan(treatedIdx(:,1)),1);
                        treatedIdx(idx) = currIndex;
                        
                        % remove it from the list
                        listCorrPx(inds==currIndex) = [];
                        inds(inds==currIndex) = [];

                        %update the currentlist to add the new data
                        currIndex = currList(1);
                        currList(1,:) = [];
                        %Quick bug fix is to check if currIndex is in inds
                        %and if not we just delete it?
                        
                        %get list of element correlated to the current element
                        list2Add  = listCorrPx{inds==currIndex}; 
                        %remove already treated cases and cases already in the list
                        list2Add(or(ismember(list2Add,treatedIdx,'row'),ismember(list2Add,currList,'row'))) = []; 

                        %find pixel that are already inside the cluster to check
                        %for correlation
                        currCluster = find(corrMask==group);
                        
                        %take a random subset of 10 pixel from the
                        %current cluster
                        try
                            r = randperm(length(currCluster,10));
                            currCluster = currCluster(r);
                        catch

                        end
                        %check that the pixel that should be added to the
                        %list are indeed correlated to the current cluster
                        [list2Add] = Core.CorrClusterMovie.enforceClusterConsistency(list2Add,...
                            currCluster,distanceMap,indsCopy);
                        %add the new element to the list
                        currList  = [currList;list2Add];

                        if isempty(currList)
                            disp(['Found ' num2str(group) ' group(s).']);
                        end

                        count = count+1;
                    end
                    %mark the last case as treated:
                    corrMask(currIndex) = group;
                    idx = find(isnan(treatedIdx(:,1)),1);
                    treatedIdx(idx) = currIndex;

                    % remove it from the list
                    listCorrPx(inds==currIndex) = [];
                    inds(inds==currIndex) = [];

                    %if exit the first while loop, the first group is complete, we need to
                    %increment the group number as we are supposed to have treated all
                    %cases.
                    group = group+1;

                end

                count=count+1;

                if count >= safeCount
                    error('something went wrong')
                end
            end   
        end
        function [neighbor] = findNeighbor(idx,dim,r)
            %function to find neighboring pixel in a given radius of a
            %central pixel given by idx. The function also makes sure that
            %we do not have indices outside the images.
            iIdx = idx(1)-r:idx(1)+r;
            jIdx = idx(2)-r:idx(2)+r;

            iIdx = iIdx(iIdx>=1 & iIdx<=dim(1));
            jIdx = jIdx(jIdx>=1 & jIdx<=dim(2));
           
            neighbor = combvec(iIdx,jIdx)';    
        end
        
        function [ind2Add] = enforceClusterConsistency(ind2Add,currCluster,distanceMap,indsCopy)
            %function to enforce consistency within the cluster by checking that
            %to-be-added element are correlated with a random subset of the cluster
            nInd = length(ind2Add);
            nClust = length(currCluster);
            newList = ind2Add;
            if ~isempty(currCluster)
                %reshape input to compare correlation
                id2Add  = find(ismember(indsCopy,newList));
                newList = repmat(id2Add,1,length(currCluster))';
                newList = newList(:);
                currCluster = find(ismember(indsCopy,currCluster));
                currCluster = repmat(currCluster,nInd,1);

                %combine input to extract value from distance map
                subs = [newList,currCluster];    
                inds = sub2ind(size(distanceMap),subs(:,1),subs(:,2));
                corrVal = distanceMap(inds);
                %reshape corrVal to test for each cluster
                corrVal = reshape(corrVal,nClust,nInd);
                %we are more relaxed since pixel far apart might be tested
                corrTest = corrVal<0.5;
                %we want new data point to be correlated to 80% of the subset from the
                %existing group
                idx2delete = sum(corrTest,1)./size(corrTest,1)<0.7;
   
                %delete the data not matching the requested threshold
                ind2Add(idx2delete) = [];    
                %we dont want any duplicate of index (although it should
                %not happen)
                ind2Add = unique(ind2Add);

            end
        end
    end
    methods(Access = 'private')
        function [run] = checkDrift(obj)
            
            path = [obj.raw.movInfo.Path filesep 'DriftCorr' filesep 'corrData.mat'];
            %test if drif corrected data exist
            test = isfile(path);
            %if drift exist we dont want to run again and vice versa
            run  = ~test;
            
            
            
        end
        
        
        
        
    end
        
end
