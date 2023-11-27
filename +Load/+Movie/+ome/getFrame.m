function [ data ] = getFrame( frameInfo, movieInfo, frames )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

assert(isvector(frames),'Frames must be a vector');
assert(isstruct(frameInfo),'info must be a structure');
assert(max(frames)<= min(movieInfo.maxFrame), 'your request exceeds the total number of frames');
    

    p2file = [movieInfo.Path filesep frameInfo.File];
    tObj   = Tiff(p2file,'r');
   
    tObj.setDirectory(1)
    data = zeros([movieInfo.Width, movieInfo.Length, length(frames)]);
    for i = 1:length(frames)
        currentDir = tObj.currentDirectory();
        data(:,:,i) = tObj.read;
        tObj.setDirectory(currentDir + 1);
        
    end
        
    
end
    
function [mov,oldPath] = loadSingleCam(movieInfo,frameInfo,frames,camIdx)
    ImL      = movieInfo.Length; % Length of the frame
    ImW      = movieInfo.Width; %width of the Frame
    nFramres = length(frames); %frame is frame2load
    mov = uint16(zeros(ImL,ImW,nFramres));
    oldPath = '';
    path2omes = movieInfo.Path;
    h = waitbar(0,'Loading a camera');
    steps = nFramres;
    for fi = 1:nFramres
        i      = camIdx(fi);
        f2load = frameInfo(i).File;
        p2file = [path2omes filesep f2load];

        if ~strcmp(oldPath,p2file)

            if exist( 'tObj', 'var' )
                tObj.close
            end
            tObj   = Tiff(p2file,'r');
            oldPath = p2file;

        end

        TifDir = str2double( frameInfo(i).IFD ) + 1;
        %     disp(TifDir)
        tObj.setDirectory(TifDir)
        mov(:,:,fi) = tObj.read;
        waitbar(fi / steps)
    end
    close (h);
    tObj.close;
end
