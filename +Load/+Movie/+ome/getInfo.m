function [frameInfo, movieInfo ] = getInfo( path2file )
%GETINFO receives as input the path to the file and gives back all
%infromation about the frames, the movie and total number of frames
%   Detailed explanation goes here
warning('off','all')
tObj = Tiff(path2file,'r');

movieInfo.Width  = tObj.getTag(256);
movieInfo.Length = tObj.getTag(257);
[folder,file,ext] = fileparts(path2file);
movieInfo.Path   = folder;

assert(tObj.currentDirectory == 1)
header = tObj.getTag(270);

idx = strfind(header,'[');
idx1 = strfind(header,',');
movieInfo.maxFrame = str2double(header(idx(1)+1:idx1(1)-1));

frameInfo.File = [file ext];

