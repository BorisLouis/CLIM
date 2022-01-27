classdef CorrClusterMovie < Core.Movie
    %General definition of a movie (master class from which many other
    %movie type will inherit
    %The Movie hold information and path to the movie but does not store
    %any data inside. Methods allow to display and do stuff with the data.
    
    properties (SetAccess = 'private')
        deconvFunction
        pathRes
        corrRelation
        corrMask
        clustEval
        corrClustMap
        cleanMask
        sil
        traceData
        
        clustLoc
        thresholdScan
        
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
                %1 Check if drift was already corrected
                [run] = obj.checkDrift;

                if or(run,strcmpi(obj.info.runMethod,'run'))
                    correlationInfo.corrSz = 75; %in px. Radius of the ROI used for correlation
                    %correlation function
                    correlationInfo.driftPeriod = 10; %in Frame, Number of frame that are averaged
                    correlationInfo.maxDrift = 20;
                    %for driftCalculation ==> 1 mean that drift is calculated for each frame
                    scalingFactor = 1;%Used for interpolation in sub-pixel Drift correction 
                    %objects of interest
                    
                    disp('Loading Data ...')
                    %h = waitbar(0,'Loading Frames');
                    nFrame = obj.raw.movInfo.maxFrame;
                    %allFrames = zeros(obj.raw.movInfo.Length,obj.raw.movInfo.Width,nFrame);
                    
                    disp(['Loading Frames ...'])
                    allFrames = obj.getFrame;
%                     for i = 1:nFrame
% 
%                         allFrames(:,:,i) = uint16(obj.getFrame(i));
%                         waitbar(i/(nFrame),h,...
%                             ['Loading Frames ' num2str(i) '/' num2str(nFrame)]);
% 
%                     end
%                     close(h);

                    disp('====>DONE<=====')
                    disp('Correcting drift...')
                    %fix Drift
                    [corrData,Drift] = PreProcess.CorrelationDrift(double(allFrames),scalingFactor,correlationInfo);

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
                    AC = autocorr(currData);
                    medAutoCorr(i,j) = AC(2);
                    
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
            
            obj.deconvFunction.Curve = deconvFunc.Data;
            obj.deconvFunction.Smoothed = deconvFunc.smoothed;
            
            figure
            subplot(1,2,1)
            imagesc(delMask(:,:,1))
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
            [run] = obj.checkCorrRelation(size(data));
            
            if or(run,strcmpi(obj.info.runMethod,'run'))
            
                data = double(data);

                r = 1;
                neigh = 8;      
               
                %2 Get correlation relationship between pixels
                [corrRelation]  = corrAnalysis.getCorrRelation(data,r,neigh);
  
                obj.corrRelation = corrRelation;
                fileName = [obj.pathRes filesep 'corrRelation.mat'];
                save(fileName,'corrRelation');
                          
            else
                
                disp('found processed data, loading from there');
                fileName = [obj.pathRes filesep 'corrRelation.mat'];
                tmp = load(fileName);
                corrR = tmp.corrRelation;
                obj.corrRelation = corrR;
                corrRelation = corrR;
                
                disp('Loading DONE');
            end
            
        end
                
        function [corrMask] = getCorrelationMask(obj,data,corrInfo)
            %Function that get correlation mask by navigating through the
            %pixel correlation map and linking one by one pixel that are
            %correlated together hence creating a map of group of pixel
            %that are correlated
            
            %check if previous data was found
            [run] = obj.checkCorrMask(size(data),corrInfo.thresh);
            
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
            subplot(1,2,2)
            RGBIM = label2rgb(corrMask,'colorcube','k','shuffle');
            imagesc(RGBIM);
            axis image
        
        end
        
        function [allCorrMask,thresh2Use] = findOptimalThreshold(obj,data,thresh)
            
            %check if previous data was found
            [run] = obj.checkThreshScan;
            
            if or(run,strcmpi(obj.info.runMethod,'run'))
                corrInfo.r = 1; %radius for checking neighbor
                corrInfo.neighbor = 8; %4 or 8 (8 means that diagonal are taken too)
                
%                 noise2Add = randn(size(data))*10;
%                 
%                 data2Save = uint16(int16(data)+int16(noise2Add));
%                 
%                 [corrRelation] = obj.getPxCorrelation(data2Save,corrInfo);
                [corrRelation] = obj.getPxCorrelation(data,corrInfo);
                allCorrMask = zeros(size(data,1),size(data,2),length(thresh));
                threshold = zeros(length(thresh),1);
                treatedArea = zeros(length(thresh),1);

                %scan threshold and determine goodness metric based on
                %silhouette
                rearrangeData = reshape(data,[size(data,1)*size(data,2),size(data,3)]);
                for i = 1:length(thresh)
                    corrInfo.thresh = thresh(i);
                    [corrM] = obj.getCorrelationMask(data,corrInfo);

                    label = reshape(corrM,[size(corrM,1)*size(corrM,2),1]);

                    %Calculate silhouette of clusters
                    tmpRearrangeData = rearrangeData;
                    tmpRearrangeData(label==0,:) =[];
                    label(label==0) = [];
                    silDist = silhouette(double(tmpRearrangeData),label,'Correlation');
       
                    sil(i) = mean(silDist);

                    threshold(i) = corrInfo.thresh;

                    %calculate % of data treated   
                    treatedArea(i) = sum(corrM(:)>0)./numel(corrM);

                    allCorrMask(:,:,i) = corrM;

                end
                % find optimal threshold
                optMetric = sil(:).*treatedArea(:);
                optMetric(isnan(optMetric)) =0;
                figure
                hold on
                plot(threshold,optMetric)
                xlabel('Correlation threshold')
                ylabel('Mean Silhouette coefficient');
                axis square
                box on

                guess.sig = 0.3;
                guess.mu = threshold(optMetric==max(optMetric));
                [FitPar,fit] = Gauss.gauss1D(optMetric,threshold,guess);

                plot(threshold,fit);

                thresh2Use = FitPar(2);
                
                %save data
                threshScan.threshold = threshold;
                threshScan.sil = sil;
                threshScan.optMetric = optMetric;
                threshScan.optMetricType = 'Sil&treatedArea';
                threshScan.treatedArea = treatedArea;
                threshScan.corrMask = allCorrMask;
                threshScan.optThresh = thresh2Use;
                
                fileName = [obj.pathRes filesep 'thresholdScan.mat'];
                save(fileName,'threshScan');
            
            else
                fileName = [obj.pathRes filesep 'thresholdScan.mat'];
                
                tmp = load(fileName);
                obj.thresholdScan = tmp.threshScan;
                
                allCorrMask = tmp.threshScan.corrMask;
                thresh2Use = tmp.threshScan.optThresh;
                
            end
            
        end
                
        function [cleanMask,silMap] = cleanCorrMask(obj,data)
            corrM = obj.corrMask.raw;
            if isempty(obj.traceData)
                traceD = obj.getAllTraces(data);
            else
                traceD = obj.traceData;
            end
                traces = [traceD.trace];
            
            clusterCorr = corrcoef(traces);
            
            stats = regionprops(corrM,'Centroid');
           
            silMap = zeros(size(corrM));
                        
            label = corrM(corrM>0);
            %px2Move = table(zeros(size(label)),zeros(size(label)),zeros(size(label)),'VariableNames',{'currClust','Sil','bestClust'});
            lab = unique(label);
            secondBestPerCluster = zeros(size(lab));
            for i = 1:length(lab)
                %gather current cluster
                idx = lab(i);
                currClust = idx;
                [row,col] = find(corrM==currClust);
                
                %get the index of top 10 closest cluster
                currCenter = stats(idx).Centroid;
                dist = cellfun(@(x) abs(x-currCenter),{stats.Centroid},'UniformOutput',false);
                euclDist = cellfun(@(x) sqrt(sum(x.^2)),dist);
                [~,idx2Clust] = mink(euclDist,10); 
                %previously getting top 10 correlated
                %[corr,idx2Clust] = maxk(clusterCorr(currClust,:),10);
                corr = clusterCorr(currClust,idx2Clust);
                traces2Comp = traces(:,idx2Clust);
                secondBestPerPixel = zeros(size(row));
                for j = 1:length(row)
                    currentTrace   = squeeze(data(row(j),col(j),:));
                    tmpTraces = [traces2Comp currentTrace];
                    corrData = corrcoef(double(tmpTraces));

                    currentCorrData = corrData(end,:);
                    currentCorrData(currentCorrData==1) = -inf;
                    
                    
                    [maxCorr,idx2MaxCorr] = maxk(currentCorrData,2);
                    %get correlation to its own cluster (where correlated
                    %cluster gave 1 in correlation)
                    correlation2Cluster = currentCorrData(corr==1);
                    
                    
                    [secondBestCorr,id] = max(maxCorr(maxCorr ~= correlation2Cluster));
                   % secondBestClust = idx2Clust(idx2MaxCorr(maxCorr==secondBestCorr));
                    silMap(row(j),col(j)) = (correlation2Cluster - secondBestCorr)/max([correlation2Cluster,secondBestCorr]);
                    
                    secondBestPerPixel(j) = idx2Clust(idx2MaxCorr(maxCorr==secondBestCorr));

%                     px2Move(i,:).currClust = currClust;
%                     px2Move(i,:).Sil  = silMap(row(i),col(i));
%                     if silMap<-0.2
%                         px2Move(i,:).bestClust = secondBestClust;
%                     else
%                         px2Move(i,:).bestClust = currClust;
%                     end

                end
                
                x = unique(secondBestPerPixel);
                bestMatchCount = groupcounts(secondBestPerPixel);
                [~,ii] = max(bestMatchCount);
                secondBestPerCluster(i) = x(ii);

            end
           
            
            cleanMask = corrM;
            cleanMask(silMap<0.2) = 0;
            obj.cleanMask = cleanMask;
            obj.sil.map = silMap;
            obj.sil.secondBest = secondBestPerCluster;
            
            figure 
            subplot(1,2,1)
            imagesc(cleanMask)
            colormap('colorcube')
            axis image
            
            subplot(1,2,2)
            imagesc(silMap)
            colormap('jet')
            axis image
          
            
            fileName = [obj.pathRes filesep 'SilCleanedCorrMask.mat'];
            save(fileName,'cleanMask');
            
            fileName = [obj.pathRes filesep 'SilhouetteMap.mat'];
            save(fileName,'silMap');
            
        end             
        function plotContour(obj,data,corrM)
            
            switch nargin
                case 2
                    assert(~isempty(obj.corrMask),'Need to run getCorrelationMask before plotting contour');
                    corrM = obj.corrMask.raw;
                case 3
                    assert(ismatrix(corrM),'correlation mask needs to be a matrix');
                otherwise 
                    error('Too many input arguments')
            end
            
            meanIm = mean(data,3);
            
            contour = cell(length(unique(corrM)),1);
            n = 1;
            for i = unique(corrM(:))'
                currMask = corrM == i;
                if not(all(currMask(:)==0))
                    contour{n} = Plotting.bwperimtrace(corrM == i,[1 size(corrM,2)],[1 size(corrM,1)]);
                    n=n+1;
                end
            end

            figure
            hold on
            imagesc(meanIm)
            for i = 1:size(contour,1)
                              

                %contour = bwboundaries(corrMaskCopy);

                plot(contour{i}{1}(:,1),contour{i}{1}(:,2),'w','LineWidth',2)

            end
            axis ij
            
            xlim([1 size(meanIm,2)])
            ylim([1 size(meanIm,1)]);
            
        end
        
        function plotClusterTraces(obj,data,idx)
            
            corrM = obj.corrMask.raw;
            [row,col] = find(corrM==idx);
            trace = zeros(length(row),size(data,3));
            figure 
            hold on
            for i = 1:length(row)

                trace(i,:) = data(row(i),col(i),:);
                plot(trace(i,:));
            end
            
        end
        
        function [clustEval1,relNum,corrClustMap] = evalCluster(obj,corrMask,data)
            
            [clustEval1,relNum,corrClustMap] = corrAnalysis.evalClusters(corrMask,data);
            relData{1} = relNum;
            label{1}   = ['Method' '-pseudoClust'];
            corrAnalysis.compareClusters(relData,label);

            obj.clustEval = clustEval1;
            obj.corrClustMap = corrClustMap;
                 
        end
        
        function [traceD] = getAllTraces(obj,data)
            assert(~isempty(obj.corrMask),'no corrMask found, please run getCorrMask first');
            corrM  = obj.corrMask.raw;
            traces = zeros(max(corrM(:)),size(data,3));
            pos    = zeros(max(corrM(:)),2);
            nClust = max(corrM(:));
            traceD = struct('trace',0,'clustPos',0,...
                'nPx',0,'clustIdx',0,'method',0);
            traceD(nClust).trace = 0;
            
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
                
                traceD(i).trace = mean(tmpTrace,2);
                traceD(i).clustPos = id;
                traceD(i).nPx = numel(id);
                traceD(i).clustIdx = i;
                traceD(i).method = 'mean';
                                            
            end
            
            obj.traceData = traceD;
            %save output
            fileName = [obj.pathRes filesep 'traceData' num2str(obj.corrMask.corrThresh) '.mat'];
            save(fileName,'traceD');
                      
        end
        
        function [corrOutput] = generateResults(obj)
            assert(~isempty(obj.traceData),'Need to run getAllTraces Firts');
            
            
            examplaryFrame = obj.getFrame(round(obj.raw.maxFrame/2));
            results = obj.traceData;
            corrMetrics = obj.clustEval;
            mask = obj.corrMask.raw;
            silM = obj.sil.map;
            secondBestClust = obj.sil.secondBest;
            nClust = length(unique(mask(mask>0)));
            assert(length(results)==nClust,'Something is inconsistent in the data');
            
            for i = 1 : nClust
                % get indices of traces 
                cMask = mask==i;
                meanSil = mean(silM(cMask));
                results(i).meanSil = meanSil;
                results(i).meanCorr = corrMetrics.meanCorr(i);
                results(i).stdCorr = corrMetrics.std(i);
                results(i).minCorr = corrMetrics.minCorr(i);
                results(i).Intensity = corrMetrics.Intensity(i);
                results(i).Amp = corrMetrics.Amp(i);
                results(i).meanInterClustCorr = corrMetrics.meanInterClustCorr(i);
                results(i).maxInterClustCorr  = corrMetrics.maxInterClustCorr(i);
                results(i).secondBestClust = secondBestClust(i);
                
            end
            
            corrOutput.deconvFunction = obj.deconvFunction;
            corrOutput.path = obj.raw.fullPath;
            corrOutput.frames = obj.raw.maxFrame;
            corrOutput.Image  = examplaryFrame;
            corrOutput.corrMap = obj.corrRelation.corrMap;
            corrOutput.corrMask = obj.corrMask.raw;
            corrOutput.rawNCluster = obj.corrMask.rawNCluster;
            corrOutput.threshold = obj.corrMask.corrThresh;
            corrOutput.method = obj.corrMask.method;
            corrOutput.cleanMask = obj.cleanMask;
            corrOutput.ROI = obj.info.ROI;
            
            if isfield(obj.info,'ROIUsed')
                corrOutput.ROIUsed = obj.info.ROIUsed;
            end
            corrOutput.silMap    = obj.sil.map;
            corrOutput.corrClustMap = obj.corrClustMap;
            corrOutput.results   = results;
            
            fileName = [obj.raw.movInfo.Path filesep 'corrAnalysisResults-' num2str(obj.corrMask.corrThresh) '.mat'];
            save(fileName,'corrOutput');

        end
        
        function saveMovie(obj,data,fps)
            path2File = obj.pathRes; 
            filename=sprintf('%s%smovie.%s', path2File,'\','avi');
            
            v = VideoWriter(filename,'Indexed AVI');
            v.Colormap = colormap(gray);
            v.FrameRate = fps;
            open(v);
            %convert to uint8
            data = double(data);
            normData = data./max(data(:))*255;
            normData(normData<0) = 0;
            
            writeVideo(v,uint8(normData(:,:,1:10:end)))
            close(v);
            
          
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
            deconvolveFunc = smooth(meanData,0.05,'rloess');
            correctedData = data;
            for i =1:size(data,1)
                for j=1:size(data,2)
                     currentData = squeeze(data(i,j,:));
                     [~,r] = deconv(currentData,meanData);
                     cleanData = r+mean(currentData);
                     correctedData(i,j,:) = cleanData; 
                end
            end
            
            deconvFunc.Data = meanData;
            deconvFunc.smoothed = deconvolveFunc;
            
        end
        
        function [RGBIM] = getImageFromMask(Mask,color)
            
            RGBIM = label2rgb(Mask,color,'k','shuffle');

            figure
            imagesc(RGBIM)
            axis image
            
            figure
            imagesc(RGBIM)
            axis image
            
            hold on
            
            centers = regionprops(Mask,'Centroid');
            
            for i = 1:length(unique(Mask(Mask>0)))
                 x = centers(i).Centroid(1);
                 y = centers(i).Centroid(2);
                 
                 text(x,y,['.',num2str(i)],'FontSize',8,'Color','white')
                
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
        
        function [run] = checkCorrRelation(obj,dim)
           
            path = [obj.pathRes filesep 'corrRelation.mat'];
            %test if file exist
            test = isfile(path);
            
            if test
                tmp = load(path);
                dimCorrRel = size(tmp.corrRelation.corrMap);
                
                if ~all(dimCorrRel == dim(1:2))
                    test = false;
                end                
            end
                
            run = ~test;
            
        end
          
        function [run] = checkCorrMask(obj,dim,thresh)
           
            path = [obj.pathRes filesep 'corrMask' num2str(thresh) '.mat'];
            test = isfile(path);
            if test%check if the mask has the correct size
                mask = load(path);
                mask = mask.cMask.raw;
                
                if ~all(size(mask)==dim(1:2))
                    test = false;
                end
            end
            %test if file exist
           
            run = ~test;
            
        end
        
        function [run] = checkThreshScan(obj)
            path = [obj.pathRes filesep 'thresholdScan.mat'];
            %run if path does not exist
            run = ~isfile(path);
        end
        
    end
        
end
