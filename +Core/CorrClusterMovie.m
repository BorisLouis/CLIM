classdef CorrClusterMovie < Core.Movie
    %General definition of a movie (master class from which many other
    %movie type will inherit
    %The Movie hold information and path to the movie but does not store
    %any data inside. Methods allow to display and do stuff with the data.
    
    properties (SetAccess = 'private')
        listCorrPx
        indCorrPx
        corrMask
        nCluster
        method
    end
        
    methods
        function obj = CorrClusterMovie(raw,info)
            
            obj = obj@Core.Movie(raw,info);
            
            if ~isfield(info,'driftCorr')
                info.driftCorr = true;
            end
            
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
          
        function [corrData] = loadFrames(obj,frames)
            %simple method to load all requested frames and allow to take
            %roi
            assert(length(frames)>=100,'Frames requested is lower than 100 please process at least 100 fr');
            fr = obj.checkFrame(frames,obj.raw.movInfo.maxFrame);
                
     
            frame = obj.getFrame(1);
            if obj.info.ROI
                figure
                
                imagesc(frame(:,:,1))
                test = drawrectangle();
                ROI  = round(test.Position);    
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
        
        function [lCorrPx,indPx] = getPxCorrelation(obj,data,corrInfo)
            data = double(data);
            r = corrInfo.r;
            corrThreshold = corrInfo.thresh;
            h=waitbar(0,'Normalizing Data...');
            % #1 Normalize the data
            data2Cluster = data-min(data,[],3);
            data2Cluster = data2Cluster./max(data2Cluster,[],3);
              
            waitbar(0.1,h,'Getting correlation relationships');
            %#2 Get correlation relationship between pixels
            [lCorrPx]  = corrAnalysis.getCorrRelation(data2Cluster,r,corrThreshold);
            
            waitbar(0.4,h,'Preparing data');
            %#3 reshap pixel relationship and get corresponding indices        
            
            indPx    = (1:length(lCorrPx))';
            
            %#4 Clean data by keeping only pixel that have correlation
            %relation
            idx2Delete = cellfun(@isempty,lCorrPx);
            lCorrPx(idx2Delete) =[];
            indPx(idx2Delete) = [];
            
            obj.listCorrPx = lCorrPx;
            obj.indCorrPx = indPx;
            
        end        
        
        function [corrMask] = getCorrelationMask(obj,data,corrInfo)
            assert(~isempty(obj.listCorrPx),'correlation relation between pixel not found, please run getPxCorrelation first');
            assert(~isempty(obj.indCorrPx), 'correlation relation between pixel not found, please run getPxCorrelation first');
            
            corrThreshold = corrInfo.thresh;
            lCorrPx = obj.listCorrPx;
            inds    = obj.indCorrPx;
            
            %get distance map
            [distanceMap] = corrAnalysis.getDistanceMapFromPxList(inds,data);
             
            disp('========> Performing Pseudo-clustering <==========')
            %perform pseudo-clustering
            dim = size(data);
            [corrMask] = corrAnalysis.corrClustering(lCorrPx,inds,distanceMap,dim,corrThreshold);
            
            disp('========> DONE <==========')
            
            figure
            imagesc(corrMask)
            axis image
            colormap('hot')
    
            obj.corrMask = corrMask;
            obj.nCluster = max(corrMask(:));
            obj.method   = 'pseudoClust';
            
        end
        
        function [MLCorrMask] = getMLCorrelationMask(obj,data,MLOptions)
            nClust = MLOptions.clust2Test;
            GPU    = MLOptions.GPU;
            replicate = MLOptions.replicate;
            dist = MLOptions.dist;
            
            inds = obj.indCorrPx;
                 
            [distanceMap]      = corrAnalysis.getDistanceMapFromPxList(inds,data);
            
            if dist
                %get xy coordinates of the pixels
                [y,x] = ind2sub(size(data),inds);
                vec = [y x];
                %calculate the distance matrix between points
                pDistance = pdist(vec);
                pDistance = squareform(pDistance);
                
                %Normalize distance so it does not overtake the correlation
                pDistance = pDistance-min(pDistance(:));
                pDistance = pDistance./max(pDistance(:));
                
                %combine correlation and spatial information in a single
                %metric
                distanceMap = distanceMap + pDistance;
                 
            end
                         
            [clust,clustEval]  = corrAnalysis.clusterCorrelatedPixel(distanceMap,...
                'clust2Test',nClust,'GPU',GPU,'replicate',replicate);
            idx = max(clust,[],1);
            clust2Use = idx == clustEval.OptimalK;
            
            MLCorrMask = zeros(size(data(:,:,1)));
            for i = 1:length(inds)
                MLCorrMask(inds(i)) = clust(i,clust2Use);
            end
           
            figure
            imagesc(MLCorrMask)
            axis image
            colormap('jet')
            
            obj.corrMask = MLCorrMask;
            obj.nCluster = max(MLCorrMask(:));
            obj.method   = 'kmeanClust';
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
        
               
        function plotContour(obj,data)
            assert(~isempty(obj.corrMask),'Need to run getCorrelationMask before plotting contour');
            maxIm = max(data,[],3);
            corrM = obj.corrMask;
            figure
            hold on
            imagesc(maxIm)
            axis image
            colormap('hot')

            for i = 1:max(corrM(:))
                corrMaskCopy = corrM;

                corrMaskCopy(corrM~=i) = 0;

                contour = bwboundaries(corrMaskCopy);

                plot(contour{1}(:,2),contour{1}(:,1),'w','LineWidth',2)

            end
            axis ij
        end
        
        function plotTraces(obj,data,idx)
            
            corrM = obj.corrMask;
            [row,col] = find(corrM==idx);
            trace = zeros(length(row),size(data,3));
            figure 
            hold on
            for i = 1:length(row)

                trace(i,:) = data(row(i),col(i),:);
                plot(trace(i,:));
            end
            
        end
        
        function [traces,pos] = getIntensityTrace(obj,data)
            corrM  = obj.corrMask;
            traces = zeros(max(corrM(:)),size(data,3));
            pos    = zeros(max(corrM(:)),2);
            for i = 1 : max(corrM(:))
                % get indices of traces 
                [row,col] = find(corrM==i);
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
                
            end
        end
        
    end
    methods(Static)
        
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
