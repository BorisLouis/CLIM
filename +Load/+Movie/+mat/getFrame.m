%% Importing data from .SPE file into matlab array 'Im' in parts
%
% Allows for frame specific import of SPE files when the entire file may
% either be too large to handle or if only a specific part of the image
% stack is important. Can be used in a for loop to divide a large stack of
% images into multiple arrays of either the same or different sizes.
%
% Ex:
%
% Im=ImportSPEframes(Ni,Nf);
%
% Ni is the starting frame and Nf is the ending frame of image stack Im. An
% error is returned if the ending frame request is larger than the actual
% image stack. Setting Ni=1 and Nf=0 returns the entire image stack.
%


function Im = getFrame(filePath,frames)
    assert(ischar(filePath),'Error, wrong path Names');
    
    data = load(filePath);
    field = fieldnames(data);
    
    mov = data.(field{1});
    nF = frames(end);
    nI = frames(1);
    
   
    Im = double(mov(:,:,nI:nF));
end