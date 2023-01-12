function [BWco]=BW2Co(BW)
% BW2Co converts a binary image into the a list od [X,Y] coordinates.
%
% Input:
%   BW:         Binary image
%   xy:         Detemines the which axis is in order, i.e. if the for loops
%               run row by row (=1) or column by column (=2). 
%               (The defult = 2, column by column)
%
% Output: 
%   BWco:       List of coodinates [X,Y]
%
% Copyright (c) 2022 O. van Zwanenberg, University of Westminster
%--------------------------------------------------------------------------
[BWco(:,2), BWco(:,1)] = find(BW~=0);