classdef CorrClusterMovie < Core.Movie
    %General definition of a movie (master class from which many other
    %movie type will inherit
    %The Movie hold information and path to the movie but does not store
    %any data inside. Methods allow to display and do stuff with the data.
    
    properties (SetAccess = 'private')
        corrRelation
        corrMask       
        pathRes
        clustLoc
        deconvFunction
        
    end
        
    methods
        function obj = CorrClusterMovie(raw,info)
            
            obj = obj@Core.Movie(raw,info);
            
            if ~isfield(info,'driftCorr')
                info.driftCorr = true;
            end
            
            folder = [obj.raw.movInfo.Path filesep 'Results'];
            
            if ~isfolder(folder)
                mkdir(folder)
            end
            
            obj.pathRes = folder;
            
        end
        
        function correctDrift(obj)
            if obj.info.driftCorr
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
            if and(obj.info.driftCorr,~obj.checkDrift)
                
                path = [obj.raw.movInfo.Path filesep 'DriftCorr' filesep 'corrData.mat'];
                assert(isfile(path),'No drift corrected data found, please run correctDrift first');
                
                dataObj = matfile(path);
                [data] = dataObj.corrData(:,:,idx);
                
            else
                
                [data] = getFrame@Core.Movie(obj,idx);
                
            end
            
        end       
          
        function [corrData] = loadFrames(obj,frames,ROI)
            %simple method to load all requested frames and allow to take
            %roi
            assert(length(frames)>=100,'Frames requested is lower than 100 please process at least 100 fr');
            fr = obj.checkFrame(frames,obj.raw.movInfo.maxFrame);
            
            switch nargin
                case 2
                    ROI = [];
                case 3
                otherwise
                    error('not enought input argument')
            end
     
            frame = obj.getFrame(1);
            if and(obj.info.ROI, isempty(ROI))  
                figure
                
                imagesc(frame(:,:,1))
                test = drawrectangle();
                ROI  = round(test.Position);    
            elseif and(obj.info.ROI,~isempty(ROI))
                assert(length(ROI)==4,'ROI is expected to be 4 elements');
                
            else
                ROI = [];
                
            end
            
            frames = obj.getFrame(fr);
                    
            if isempty(ROI)
               ROI = [1 1 size(frames,2),size(frames,1)];
            end
            
            %Fix ROI
            if mod(ROI(3),2) ~=0
                ROI(3) = ROI(3)-1;
            end
            
            if mod(ROI(4),2) ~=0
                ROI(4) = ROI(4)-1;
            end
            %apply ROI
            corrData = frames(ROI(2):ROI(2)+ROI(4)-1,ROI(1):ROI(1) + ROI(3)-1,:);
       
        end
        
        function [correctedData] = deconvolve(obj,data)
            
            disp('Detecting background based on autocorrelation')
            %remove background
            medAutoCorr = zeros(size(data,1),size(data,2));
            for i = 1:size(data,1)
                for j = 1:size(data,2)
                
                    currData = double(squeeze(data(i,j,:)));
                    medAutoCorr(i,j) = median(autocorr(currData));
                    
                end
            end
            
            deleteMask = medAutoCorr>0.1;
            %clean up the binary image mask
            SE = strel('disk',4);
            deleteMask = imopen(deleteMask,SE);
            deleteMask = imfill(deleteMask,'holes');
            %keep only the largest area
            d = regionprops(deleteMask,'Area','PixelIdxList');
            [~,idx] = max([d.Area]);
            delMask = zeros(size(deleteMask));
            delMask(d(idx).PixelIdxList) = 1;
            %store the mask            
            obj.deconvFunction.Mask = delMask;
            %repeeat it in z and multiply it
            delMask = repmat(delMask,1,1,size(data,3));
            
            bkgCorrData = double(delMask).*double(data);
                     
            %deconvolve data
            [correctedData,deconvFunc] = Core.CorrClusterMovie.deconvolveFromMean(bkgCorrData);
            
            obj.deconvFunction.Curve = deconvFunc;
            
            figure
            subplot(1,2,1)
            imagesc(deleteMask(:,:,1))
            title('Background removal mask')
            axis image
            subplot(1,2,2)
            plot(deconvFunc.Data)
            hold on
            plot(deconvFunc.smoothed,'r','Linewidth',2);
            legend({'Data','Smoothed'})
            axis square
            title('Deconvolution curve')
            box on
            
        end
        

        
        function [corrRelation] = getPxCorrelation(obj,data,corrInfo)
            %This function scan the image to find pixel that are correlated
            %to other pixel that are close in space (typically check only
            %neighboring pixels
            
            %check if previous data was found
            [run] = obj.checkPxData;
            
            if or(run,strcmpi(obj.info.runMethod,'run'))
            
                data = double(data);
                r = corrInfo.r;
                h=waitbar(0,'Normalizing Data...');
%                 % #1 Normalize the data
%                 Does not affect the data so deleted this part
%                 data2Cluster = data-min(data,[],3);
%                 data2Cluster = data2Cluster./max(data2Cluster,[],3);
            

                waitbar(0.1,h,'Getting correlation relationships');
                %#2 Get correlation relationship between pixels
                [corrRelation]  = corrAnalysis.getCorrRelation(data,r);

                waitbar(0.8,h,'Saving data');
  
                obj.corrRelation = corrRelation;
                fileName = [obj.pathRes filesep 'corrRelation.mat'];
                save(fileName,'corrRelation');
                waitbar(1,h,'Finished');
                pause(1);
                close(h);
                
            else
                
                disp('found processed data, loading from there');
                fileName = [obj.pathRes filesep 'corrRelation.mat'];
                tmp = load(fileName);
                corrR = tmp.corrR;
                obj.corrRelation = corrR;
                
                disp('Loading DONE');
            end
            
        end
                
        function [corrMask] = getCorrelationMask(obj,data,corrInfo)
            %Function that get correlation mask by navigating through the
            %pixel correlation map and linking one by one pixel that are
            %correlated together hence creating a map of group of pixel
            %that are correlated
            
            %check if previous data was found
            [run] = obj.checkCorrMask(corrInfo.thresh);
            
            if or(run,strcmpi(obj.info.runMethod,'run'))
                assert(~isempty(obj.corrRelation),'correlation relation between pixel not found, please run getPxCorrelation first');
                assert(~isempty(obj.corrRelation.listPx),'correlation relation between pixel not found, please run getPxCorrelation first');
                assert(~isempty(obj.corrRelation.indPx), 'correlation relation between pixel not found, please run getPxCorrelation first');

                corrThreshold = corrInfo.thresh;
                listPx = obj.corrRelation.listPx;
                inds    = obj.corrRelation.indPx;
                meanPx   = obj.corrRelation.meanPx;
                listVal = obj.corrRelation.listVal;
                %get distance map
                %[distanceMap] = corrAnalysis.getDistanceMapFromPxList(inds,data);

                disp('========> Performing Pseudo-clustering <==========')
                %perform pseudo-clustering
                [corrMask] = corrAnalysis.corrClustering(listPx,listVal,meanPx,inds,data,corrThreshold);                
                
                %save Data 
                cMask.raw = corrMask;
                cMask.rawNCluster = max(corrMask(:));
                cMask.method = 'pseudoClust';
                cMask.corrThresh = corrThreshold;
                
                obj.corrMask = cMask;
                
                fileName = [obj.pathRes filesep 'corrMask' num2str(obj.corrMask.corrThresh) '.mat'];
                save(fileName,'cMask');

            else
                disp('Found CorrMask from previous analysis, loading from there');
                
                fileName = [obj.pathRes filesep 'corrMask' num2str(corrInfo.thresh) '.mat'];
                tmp = load(fileName);
                
                obj.corrMask = tmp.cMask;
                
                corrMask = obj.corrMask.raw;       
            end
            disp('========> DONE <==========')
            
            figure
            subplot(1,2,1)
            imagesc(corrMask)
            axis image
            colormap('jet')
        
        end
                
        function [meanTraces] = checkMask(obj,data,clust)
            mask = obj.corrMask;
            data = double(data);
            %get the requested cluster
            subMask = mask ==clust;
            subMask = bwareaopen(subMask,8);
            subMask = imfill(subMask,8,'holes');
            
            labSubMask = bwlabel(subMask);
            
            idx = max(labSubMask(:));
            
            disp(['Found ' num2str(idx) ' spatially speparated clusters with the number ' num2str(clust)]);
            
            meanTraces = zeros(idx,size(data,3));
            for i = 1:idx
                % get indices of traces 
                [row,col] = find(labSubMask==i);
                r = repmat(row,1,size(data,3));
                r = reshape(r',size(data,3)*length(row),1);
                c = repmat(col,1,size(data,3));
                c = reshape(c',size(data,3)*length(col),1);
                
                f = repmat((1:size(data,3))',length(col),1);
                idx = sub2ind(size(data),r,c,f);
                
                %get the data
                tmpTrace = data(idx);
                tmpTrace = reshape(tmpTrace,size(data,3),length(row));
                
                %check correlation
                tmpCorr = corrcoef(double(tmpTrace));
                                              
                disp(['The min correlation in the cluster is ' num2str(min(tmpCorr(:)))]);
                disp(['The max correlation in the cluster is ' num2str(max(tmpCorr(tmpCorr<1)))]);
                disp(['The median correlation in the cluster is ' num2str(median(tmpCorr(tmpCorr<1)))])
                
                meanTraces(i,:) = mean(tmpTrace,2);
                
                id = randperm(size(tmpTrace,2),5);
                
                trace2Show = tmpTrace(:,id);
                figure
                hold on
                for j = 1 : length(id)
                   
                    plot(trace2Show(:,j)./mean(trace2Show(:,j),1));
                    
                end
                plot(meanTraces(i,:)./mean(meanTraces(i,:),2),'linewidth',2,'color',[0.4 0.4 0.4]);
                xlabel('Frames')
                ylabel('Norm. Intensity')
                axis square
                box on
                title(['Cluster ' num2str(clust) '.' num2str(i) ' examplary traces and avg']);
                
            end
            
            meanCorr = corrcoef(meanTraces');
            uniCorr = unique(meanCorr);
            disp(['The spatially separated clusters are correlated together with: ' num2str(uniCorr(uniCorr<1)')]);
            
            figure
            hold on
            for i = 1:size(meanTraces,1)
                plot(meanTraces(i,:)./mean(meanTraces(i,:)))
            end
            xlabel('Frames')
            ylabel('Norm. Intensity')
            title('comp. cluster with same number')
            
            
        end
              
        function plotContour(obj,data,corrM)
            
            switch nargin
                case 2
                    assert(~isempty(obj.corrMask),'Need to run getCorrelationMask before plotting contour');
                    corrM = obj.corrMask.raw;
                case 3
                otherwise 
                    error('Too many input arguments')
            end
            
            maxIm = max(data,[],3);
            
            figure
            hold on
            imagesc(maxIm)
            axis image
            colormap('jet')

            for i = 1:max(corrM(:))
                corrMaskCopy = corrM;

                corrMaskCopy(corrM~=i) = 0;

                contour = bwboundaries(corrMaskCopy);

                plot(contour{1}(:,2),contour{1}(:,1),'w','LineWidth',2)

            end
            axis ij
        end
        
        function plotClusterTraces(obj,data,idx)
            
            corrM = obj.corrMask.clean;
            [row,col] = find(corrM==idx);
            trace = zeros(length(row),size(data,3));
            figure 
            hold on
            for i = 1:length(row)

                trace(i,:) = data(row(i),col(i),:);
                plot(trace(i,:));
            end
            
        end
        
        function [traceData] = getAllTraces(obj,data)
            assert(~isempty(obj.corrMask),'no corrMask found, please run getCorrMask first');
            corrM  = obj.corrMask.clean;
            traces = zeros(max(corrM(:)),size(data,3));
            pos    = zeros(max(corrM(:)),2);
            nClust = max(corrM(:));
            traceData = struct('trace',0,'clustPos',0,...
                'nPx',0,'method',0);
            traceData(nClust).trace = 0;
            
            for i = 1 : nClust
                % get indices of traces 
                [row,col] = find(corrM==i);
                id = find(corrM==i);
                r = repmat(row,1,size(data,3));
                r = reshape(r',size(data,3)*length(row),1);
                c = repmat(col,1,size(data,3));
                c = reshape(c',size(data,3)*length(col),1);
                
                f = repmat((1:size(data,3))',length(col),1);
                idx = sub2ind(size(data),r,c,f);
                
                tmpTrace = data(idx);
                tmpTrace = reshape(tmpTrace,size(data,3),length(row));
                
                traces(i,:) = mean(tmpTrace,2);
                pos(i,:)    = [mean(row), mean(col)];
                
                traceData(i).trace = mean(tmpTrace,2);
                traceData(i).clustPos = id;
                traceData(i).nPx = numel(id);
                traceData(i).method = 'mean';
                                            
            end
            
            %save output
            fileName = [obj.pathRes filesep 'traceData' num2str(obj.corrMask.corrThresh) '.mat'];
            save(fileName,'traceData');
                      
        end
        
        function [RGBIM] = getImageFromMask(obj,Mask,color)
            
            colors = colormap(color);
            colors = [0,0,0; colors];

            RGBIM = zeros(size(Mask,1),size(Mask,2),3);
            
            [row,col] = find(Mask==0);
            r = repmat(row,1,size(RGBIM,3));
            r = reshape(r',size(RGBIM,3)*length(row),1);
            c = repmat(col,1,size(RGBIM,3));
            c = reshape(c',size(RGBIM,3)*length(col),1);

            f = repmat((1:size(RGBIM,3))',length(col),1);
            idx = sub2ind(size(RGBIM),r,c,f);
            %color background
            
            RGBIM(idx) = repmat(colors(1,:)',length(row),1);
            %for the 'real' data use everything above 25% of the range so
            %the background is more constrasted.
            idx2Use = round(0.25*length(colors));
            cropColors = colors(idx2Use:end,:);
            for i = 1:max(Mask(:))
                [row,col] = find(Mask==i);
                r = repmat(row,1,size(RGBIM,3));
                r = reshape(r',size(RGBIM,3)*length(row),1);
                c = repmat(col,1,size(RGBIM,3));
                c = reshape(c',size(RGBIM,3)*length(col),1);

                f = repmat((1:size(RGBIM,3))',length(col),1);
                idx = sub2ind(size(RGBIM),r,c,f);
                
                idx2Use = randi(size(cropColors,1));
                
                RGBIM(idx) = repmat(cropColors(idx2Use,:)',length(row),1);
                
                
            end
            
            figure
            imagesc(RGBIM)
            axis image
            
            
            
        end
        
        function [loc] = getClusterLocalization(obj,data,idx)
            cMask = obj.corrMask;
            mask = cMask.raw;
            
            images = regionprops(mask,'Image','BoundingBox');
            
            bBox = round(images(idx).BoundingBox);
            currentImage = images(idx).Image;
            
            data2Fit = double(data(bBox(2):bBox(2)+bBox(4)-1,bBox(1):bBox(1)+bBox(3)-1,:));
            
            [X,Y]= meshgrid(bBox(1):bBox(1)+bBox(3)-1,bBox(2):bBox(2)+bBox(4)-1);
            domain(:,:,1) = X;
            domain(:,:,2) = Y;
            
            loc = struct('x',zeros(size(data2Fit,3),1),'y',zeros(size(data2Fit,3),1),...
                'int',zeros(size(data2Fit,3),1),'angle',zeros(size(data2Fit,3),1));
            
            for i = 1:size(data2Fit,3)
                currentFrame = data2Fit(:,:,i);
                currentFrame(~currentImage) = 0;
                
                %fit the region with a Gaussian
                [gPar, fit] = Gauss.Gauss2D_Fit(currentFrame,domain);
%                 figure
%                 surf(currentFrame)
%                 hold on 
%                 surf(fit)
                
                loc.x(i) = gPar(5);
                loc.y(i) = gPar(6);
                loc.int(i) = gPar(1);
                loc.angle(i) = gPar(7);
                
            end 
        end
        
        function [allLoc] = getAllClusterLocalization(obj,data)
           nCluster = obj.corrMask.rawNCluster;
           allLoc = cell(nCluster,1);
           for i = 1 :nCluster
               
               loc = obj.getClusterLocalization(data,i);
               allLoc{i} = loc;
               
           end
           obj.clustLoc = allLoc;
        end
        
        function showClusterLoc(obj,idx)
            assert(~isempty(obj.clustLoc),'Please run getAllClusterLocalization');
            
            locData = obj.clustLoc{idx};
            
            pos = norm([locData.x,locData.y]);
            int = locData.int;
            angle = locData.angle;
            
            figure
            subplot(1,3,1)
            scatter(pos,int,10,'filled')
            axis square
            box on
            title('Intensity vs Position')
            
            subplot(1,3,2)
            scatter(pos,angle,10,'filled')
            axis square
            box on
            title('Position vs angle')
            
            subplot(1,3,3)
            scatter(int,angle,10,'filled')
            axis square
            box on
            title('Intensity vs Angle')
            
        end
        
    end
    methods(Static)
        function [correctedData,deconvFunc] = deconvolveFromMean(data)
            
            %calculate the mean intensity trace, ignoring zeros from
            %background deletion step
            meanData = squeeze(sum(sum(data,1),2)./sum(sum(data(:,:,1)~=0,1),2));
            deconvolveFunc = smooth(meanData,0.1,'rloess');
            correctedData = data;
            for i =1:size(data,1)
                for j=1:size(data,2)
                     currentData = squeeze(data(i,j,:));
                     [~,r] = deconv(currentData,deconvolveFunc);
                     cleanData = r+mean(currentData);
                     correctedData(i,j,:) = cleanData; 
                end
            end
            
            deconvFunc.Data = meanData;
            deconvFunc.smoothed = deconvolveFunc;
            
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
        
        function [run] = checkPxData(obj)
           
            path = [obj.pathRes filesep 'pxData.mat'];
            %test if file exist
            test = isfile(path);
            
            run = ~test;
            
        end
        
        
        function [run] = checkCorrMask(obj,thresh)
           
            path = [obj.pathRes filesep 'corrMask' num2str(thresh) '.mat'];
            %test if file exist
            test = isfile(path);
            run = ~test;
            
        end
        
    end
        
end
