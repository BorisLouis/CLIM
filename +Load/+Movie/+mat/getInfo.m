function [frameInfo,movInfo] = getInfo(path2File)
    %TODO: Extract exposure time 
    [path,file,ext] = fileparts(path2File);
    
    data = load(path2File);
    field = fieldnames(data);
    
    mov = data.(field{1});
    
    xDim = size(mov,2);
    yDim = size(mov,1);
    nFrames = size(mov,3);
      
    %get the number of frame
    maxFrame = nFrames;
    %get the exposure time
    expT     = NaN; %in ms
    %store a few thing for the rest of the code
    isMultiImage = false;
    
    isZStack = false;
    
    Cam  = 0;
    
    X = xDim;
    Y = yDim;
        
    %Store info for output
    frameInfo.File = [file ext];
    
    movInfo.Width  = X;
    movInfo.Length = Y;
    movInfo.Path   = path;
    movInfo.isMultiImage = isMultiImage;
    movInfo.isZStack = isZStack;
    movInfo.Cam = Cam;
    movInfo.expT = expT;
    movInfo.maxFrame = maxFrame;

end
