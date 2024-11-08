classdef CorrClusterMovie < Core.Movie
    %General definition of a movie (master class from which many other
    %movie type will inherit
    %The Movie hold information and path to the movie but does not store
    %any data inside. Methods allow to display and do stuff with the data.
    
    properties (SetAccess = 'private')
        deconvFunction
        pathRes
        corrRelation
        
        clustEval
        corrClustMap
        cleanMask
        sil
        traceData
        
        clustLoc
        thresholdScan
        
    end
    properties (SetAccess = 'public')
        corrMask
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
        
        function [correctedData] = deconvolve(obj,data,decThresh,doDeconvolve)
            
            disp('Detecting background based on autocorrelation')
            %basic segmentation
            d = mean(data,3);
            bw = imbinarize(d./max(d(:)),0.1);
            
            area = regionprops(bw,'area','pixelidxlist');
            %get largest area
            [~,id] = max([area.Area]);
            
            pixels = area(id).PixelIdxList;
            mask = zeros(size(d));
            mask(pixels) = 1;
            mask = imfill(mask,'holes');
            bkg  = double(data.*uint16(~bw));
            bkg(bkg==0) = NaN;
            background = squeeze(nanmedian(nanmedian(bkg,1),2));
            signal = data.*uint16(mask);
     %                      
            %deconvolve data
            if doDeconvolve
                [correctedData,deconvFunc] = Core.CorrClusterMovie.deconvolveFromMean(signal);
            else
                deconvFunc.Data = ones(1:size(data,3));
                deconvFunc.smoothed = deconvFunc.Data;
                correctedData = signal;
            end
            
            %store the mask
            obj.deconvFunction.deconv = doDeconvolve;
            obj.deconvFunction.Mask = mask;
            obj.deconvFunction.Curve = deconvFunc.Data;
            obj.deconvFunction.Smoothed = deconvFunc.smoothed;
            obj.deconvFunction.background = background;
            obj.deconvFunction.background2Sub = mean(background);
            
   
            figure
            subplot(1,2,1)
            imagesc(mask(:,:,1))
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
               
        function [corrRelation] = getPxCorrelation(obj,data)
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
                switch lower(obj.info.thresholdMode)
                    case 'none'
                        [corrMask] = corrAnalysis.corrClusteringNoThresh(obj.corrRelation,data,corrThreshold,obj.info.doPlot);                
                        cMask.method = 'noThreshold';
                    otherwise
                    %perform pseudo-clustering
                    [corrMask] = corrAnalysis.corrClustering(listPx,listVal,meanPx,inds,data,corrThreshold,obj.info.doPlot);                
                    cMask.method = 'pseudoClust';                 
                end
                
                %save Data 
                cMask.raw = corrMask;
                cMask.rawNCluster = max(corrMask(:));
                
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

                [~] = obj.getPxCorrelation(data);
                allCorrMask = zeros(size(data,1),size(data,2),length(thresh));
                threshold = zeros(length(thresh),1);
                treatedArea = zeros(length(thresh),1);

                %scan threshold and determine goodness metric based on
                %silhouette
                rearrangeData = reshape(data,[size(data,1)*size(data,2),size(data,3)]);
                Sil = zeros(length(thresh),1);
                SilClust = Sil;
                for i = 1:length(thresh)
                    corrInfo.thresh = thresh(i);
                    [corrM] = obj.getCorrelationMask(data,corrInfo);

                    label = reshape(corrM,[size(corrM,1)*size(corrM,2),1]);

                    %Calculate silhouette of clusters
                    tmpRearrangeData = rearrangeData;
                    tmpRearrangeData(label==0,:) =[];
                    label(label==0) = [];
                    %old silhouette calculation
%                     silDist = silhouette(double(tmpRearrangeData),label,'Correlation');
%        
%                     Sil(i) = mean(silDist);
                    
                    % get silhouette from clusters
                    [~,silMap] = obj.cleanCorrMask(data);
                    SilClust(i) = mean(silMap(:));
                    
                    % end addition
                    
                    threshold(i) = corrInfo.thresh;

                    %calculate % of data treated   
                    treatedArea(i) = sum(corrM(:)>0)./numel(corrM);

                    allCorrMask(:,:,i) = corrM;

                end
                % find optimal threshold
                optMetric = SilClust(:).*treatedArea(:);
                optMetric(isnan(optMetric)) =0;
                figure
                hold on
                plot(threshold,optMetric)
                xlabel('Correlation threshold')
                ylabel('Mean Silhouette coefficient');
                axis square
                box on

%                 guess.sig = 0.3;
%                 guess.mu = 0.6;
%                 [FitPar,fit] = Gauss.gauss1D(optMetric,threshold,guess);
% 
%                 plot(threshold,fit);
                [~,idx] = max(optMetric);
                thresh2Use = threshold(idx);
                
                %save data
                threshScan.threshold = threshold;
                threshScan.sil = Sil;
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
            
            if length(obj.traceData) ~= max(corrM(:))
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
            if length(lab) >1
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
               
                    corr = clusterCorr(currClust,idx2Clust);
                    traces2Comp = traces(:,idx2Clust);
                    secondBestPerPixel = zeros(size(row));
                    for j = 1:length(row)
                        currentTrace   = squeeze(data(row(j),col(j),:));
                        tmpTraces = [traces2Comp currentTrace];
                        corrData = corrcoef(double(tmpTraces));

                        currentCorrData = corrData(end,:);
                        correlation2Cluster = currentCorrData(corr==1);
                        currentCorrData(currentCorrData==1) = -inf;


                        [maxCorr,idx2MaxCorr] = maxk(currentCorrData,2);
                        %get correlation to its own cluster (where correlated
                        %cluster gave 1 in correlation)



                        [secondBestCorr,id] = max(maxCorr(maxCorr ~= correlation2Cluster));
                       % secondBestClust = idx2Clust(idx2MaxCorr(maxCorr==secondBestCorr));
                       if all(isnan(currentCorrData))
                           silMap(row(j),col(j)) = 0;
                           secondBestPerPixel(j) = NaN;
                       else
                        silMap(row(j),col(j)) = (correlation2Cluster - secondBestCorr)/max([correlation2Cluster,secondBestCorr]);

                        secondBestPerPixel(j) = idx2Clust(idx2MaxCorr(maxCorr==secondBestCorr));
                       end
  
                    end

                    x = unique(secondBestPerPixel);
                    bestMatchCount = groupcounts(secondBestPerPixel);
                    [~,ii] = max(bestMatchCount);
                    secondBestPerCluster(i) = x(ii);

                end
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
            for i = 1:10:length(row)

                trace(i,:) = data(row(i),col(i),:);
                normTrace = (trace(i,:) - min(trace(i,:)))./(max(trace(i,:)) - min(trace(i,:)));
                plot(trace(i,:));
               
            end
            disp('stop');
        end
        
        function [clustEval1,relNum,corrClustMap] = evalCluster(obj,corrMask,data)
            
            [clustEval1,relNum,corrClustMap] = corrAnalysis.evalClusters(corrMask,data);
            relData{1} = relNum;
            label{1}   = ['Method' '-pseudoClust'];
            corrAnalysis.compareClusters(relData,label);

            obj.clustEval = clustEval1;
            obj.corrClustMap = corrClustMap;
                 
        end
        
        function [traceD] = getAllTraces(obj,data,corrData,method)
            
            switch nargin
                case 2
                    method = 'mean';
                    corrData = 0;
                case 3
                    method = 'mean';
                case 4 
                otherwise
                    error('Too many input argument')
            end
            
            if ~isempty(obj.sil)
                silVal = obj.sil.map;
            else
                method = 'mean';
                warning('Intensity extraction method chosen could not be achieve because the silhouette map is missing, calculating via mean...');
                
            end
            
            assert(~isempty(obj.corrMask),'no corrMask found, please run getCorrMask first');
            corrM  = obj.corrMask.raw;
            traces = zeros(max(corrM(:)),size(data,3));
            pos    = zeros(max(corrM(:)),2);
            nClust = max(corrM(:));
            traceD = struct('trace',0,'corrTrace',0,'clustPos',0,...
                    'nPx',0,'clustIdx',0,'method',0);
                
            if nClust ~= 0
                
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
                    
                    if length(corrData)~=1
                        tmpCorrTrace = corrData(idx);
                        tmpCorrTrace = reshape(tmpCorrTrace,size(data,3),length(row));
                    else
                        tmpCorrTrace = 0;
                    end


                    pos(i,:)    = [mean(row), mean(col)];

                    switch lower(method)
                        case 'mean'
                             traceD(i).trace = mean(tmpTrace,2);
                             traceD(i).corrTrace = mean(tmpCorrTrace,2);
                        case 'silweigth'
                             silWeight = silVal(id);
                             tmpTrace(:,silWeight<0.05) = [];
                             if length(corrData)~=1
                                tmpCorrTrace(:,silWeight<0.05) = [];
                             end
                             silWeight(silWeight<0.05) = [];
                             silWeight = repmat(silWeight',size(tmpTrace,1),1);
                             traceD(i).trace = sum(double(tmpTrace).*silWeight,2)./(sum(silWeight(1,:)));
                             if length(corrData)~=1
                                traceD(i).corrTrace = sum(double(tmpCorrTrace).*silWeight,2)./(sum(silWeight(1,:)));
                             end
                        case 'silthresh'
                            silWeight = silVal(id);
                            traceD(i).trace = mean(tmpTrace(:,silWeight>0.2),2);
                            if length(corrData)~=1
                                traceD(i).corrTrace= mean(tmpCorrTrace(:,silWeight>0.2),2);
                            end
                    end

                    traceD(i).clustPos = id;
                    traceD(i).nPx = numel(id);
                    traceD(i).clustIdx = i;
                    traceD(i).method = method;

                end
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
            else
                corrOutput.ROIUsed = [];
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
            %sum in 2D across time (spatial averaged time trace) and divide
            %by nnz (number of nonzero element) in one of the slice to
            %avoid weighing in the zeros from backgound
            meanData = squeeze(sum(sum(data,1),2)./nnz(data(:,:,1)));
            deconvolveFunc = smooth(meanData,0.05,'rloess');
            correctedData = data;
            for i =1:size(data,1)
                for j=1:size(data,2)
                     currentData = squeeze(data(i,j,:));
                     if ~all(currentData==0)
                         [~,r] = deconv(currentData,meanData);
                         cleanData = r+mean(currentData);
                         correctedData(i,j,:) = cleanData;
                     else
                         correctedData(i,j,:) = currentData;
                     end
                         
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
