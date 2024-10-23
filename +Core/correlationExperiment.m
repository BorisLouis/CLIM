classdef correlationExperiment < handle
    
    properties
        path
        ext
        corrMovies
        info
        corrMasks
        nClusters
    end
    
    methods
        
        function obj = correlationExperiment(folder2Data,info)
            
            %correlationExperiment Construct an instance of this class
            %   Detailed explanation goes here
            obj.path = folder2Data.path;
            obj.ext  = folder2Data.ext;
            
            %check info present in info to make sure everything is provided
            if ~isfield(info,'driftCorr')
                info.driftCorr = true;
                warning('No driftcorr info provided, correcting drift by default');
            end
            
            if ~isfield(info,'corrInfo')
                info.corrInfo.thresh = 0.5;
                warning('No correlation information provided, set r=2 and thresh= 0.3 as default value');
                
            end
            
            if ~isfield(info,'frame2Process')
                info.frame2Process = 1:1000;
                warning('No frame2Process information provided, set 1:1000 as default')
            end
          
            obj.info = info;
        end
        
        %Set function
        function set.path(obj, path)
            assert(ischar(path), 'Path should be given as a string');
            assert(isfolder(path), 'The path given is not a folder, ZCalibration expect a folder. In the folder it is expected to find separate folder for each zCalMovie.')
            
            obj.path = path;
            
        end
        
        
        function retrieveMovies(obj)
            
            folder2Mov = dir(obj.path);
            folder2Mov = folder2Mov(cell2mat({folder2Mov.isdir}));
            %loop through the content of the directory
            count = 0;
            for i = 3:size(folder2Mov,1)
                %Check if the directory
                folderPath = [folder2Mov(i).folder filesep folder2Mov(i).name];
                file2Analyze = Core.Movie.getFileInPath(folderPath,obj.ext);

                if ~isempty(file2Analyze)
                    count = count+1;
                    file.path = file2Analyze.folder;
                    file.ext  = obj.ext;
                    tmp = Core.CorrClusterMovie(file, obj.info);

                    if count == 1
                        tmp.giveInfo;
                    else
                        tmp.info = obj.corrMovies.(['mov' num2str(1)]).getInfo; 
                    end
                  
                    obj.corrMovies.(['mov' num2str(i-2)]) = tmp;


                else

                    warning([folder2Mov(i).folder filesep folder2Mov(i).name ' did not contain any ' obj.ext ' file and is therefore ignored']);

                end

            end
            

            if isempty(obj.corrMovies)              
               error(['No %s was found. Please check:\n',...
                   '1) that you gave the correct file extension.\n',...
                   '2) that you gave the path of a folder containing folders containing movies with the given extension'],obj.ext);        

            end
            disp('=======> DONE ! <========')
        end
        
        
        function getCorrMask(obj)
            
            fieldsN = fieldnames(obj.corrMovies);
            %Extraction of Data
            nfields = numel(fieldsN);
            corrInfo = obj.info.corrInfo;
            f2Process = obj.info.frame2Process;
            threshold = obj.info.threshold;
            for i = 1: nfields
                
                disp(['Retrieving data from fluctuating file ' num2str(i) ' / ' num2str(nfields) ' ...']);
                currentCorrMov = obj.corrMovies.(fieldsN{i});
                
                currentCorrMov.correctDrift;
                frames = Core.Movie.checkFrame(f2Process, currentCorrMov.raw.movInfo.maxFrame);
                data2Use = currentCorrMov.loadFrames(frames);
                
                %% Scanning threshold
                if strcmpi(obj.info.thresholdMode,'auto')
                    testROIRadius = 64;
                    center = [round(size(data2Use,1)/2), round(size(data2Use,2)/2)];
                    if obj.info.ROI ==false
                        currentCorrMov.info.ROIUsed = [];    
                    else
                        ROICorrData = data2Use;
                        currentCorrMov.info.ROIUsed = ROI;
                    end
                    
                    try

                        ROICorrData = data2Use(center(1)-testROIRadius:center(1)+testROIRadius-1,center(2)-testROIRadius:center(2)+testROIRadius-1,:);
                    catch except
                        if strcmp(except.identifier, 'MATLAB:badsubscript')
                            ROICorrData = data2Use;
                        end
                    end

                    thresh = obj.info.minCorr:obj.info.stepCorr:obj.info.maxCorr;

                    [allCorrMask,threshold2Use] = currentCorrMov.findOptimalThreshold(ROICorrData,thresh);
                end
                %% Calculate final corrMask
                switch lower(obj.info.thresholdMode)
                    case 'auto'
                        corrInfo.thresh = threshold2Use;
                    case 'fixed'
                        corrInfo.thresh = threshold;
                    case 'none'
                    corrInfo.thresh = obj.info.minCorr;
                end
                %get px correlation
                [corrRelation] = currentCorrMov.getPxCorrelation(data2Use);
                %get the mask using the optimal threshold
                [corrMask] = currentCorrMov.getCorrelationMask(data2Use,corrInfo);

                obj.corrMasks(i).map = corrRelation.corrMap;
                obj.corrMasks(i).mask= corrMask;
                obj.corrMasks(i).data = mean(data2Use,3);
                obj.nClusters(i) = max(corrMask(:));

            end
        end
        
        function s = saveobj(obj)
            s.path = obj.path;
            s.ext  = obj.ext;
            s.corrMovies = obj.corrMovies;
            s.info = obj.info;
            s.corrMasks = obj.corrMasks;
            s.nClusters = obj.nClusters;
            
        end
             
    end
    
    
end