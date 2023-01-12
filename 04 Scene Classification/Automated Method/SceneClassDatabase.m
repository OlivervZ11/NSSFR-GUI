% Scene Classification of the database

% O. van Zwanenberg (June 2020)
% UNIVERSITY OF WESTMINSTER 
%              - COMPUTATIONAL VISION AND IMAGING TECHNOLOGY RESEARCH GROUP

clc; close all; clear all;
load('SCNet.mat')

%--------------------------------------------------------------------------

path = uigetdir;
filePattern = fullfile(path, '*.tif');
jpegFiles = dir(filePattern);
tic
names={jpegFiles.name}';

% imnumber stores the number of files that have been read
imnumber=size(names,1);

% Large loop of functions to extract edges and MTFs from each image in imfiles structure

Indoor=0;
ManMade=0;
Nature=0;

parfor A=1:imnumber
    i = fullfile(path,names{A});
    Im=imread(i);
    Im=im2uint8(Im);
    
    image = imresize(Im,[227 227]);
    
    [category, score] = trainedNet.classify(image);

    category=char(category(1));
    
    if strcmp(category,'Indoor')
        Indoor=Indoor+1;
    elseif strcmp(category,'ManMade')
        ManMade=ManMade+1;  
    elseif strcmp(category,'Natural')
        Nature=Nature+1;
    end  
end

Indoor=(Indoor/imnumber)*100;
ManMade=(ManMade/imnumber)*100;
Nature=(Nature/imnumber)*100;

sceneclass=[Indoor,ManMade,Nature];
toc
clearvars -except sceneclass