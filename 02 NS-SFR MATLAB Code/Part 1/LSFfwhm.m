function [LSFHPW]=LSFfwhm(LSF)
% LSFhpw converts the ESF to the LSF, and then takes Full Width Half Max
%
% Input: 
%       LSF       -    The Line Spread Function of interest
% Output: 
%       LSFHPW    -    The LSF Half Peak Width
%
% Copyright (c) 2020 O. van Zwanenberg
% UNIVERSITY OF WESTMINSTER 
%              - COMPUTATIONAL VISION AND IMAGING TECHNOLOGY RESEARCH GROUP
% Director of Studies:  S. Triantaphillidou
% Supervisory Team:     R. Jenkin & A. Psarrou
 
[m,n] = size(LSF);
if m>n
    a=1;
else
    a=2;
end
x=1:size(LSF,a);

% Caculate the Half Peak Width

% Find the half height - midway between min and max LSF values
halfHeight = (min(LSF) + max(LSF)) / 2;

% Find left edge
index1 = find(abs(LSF) >= abs(halfHeight), 1, 'first');
x1 = x(index1);
% Find right edge
index2 = find(abs(LSF) >= abs(halfHeight), 1, 'last');
x2 = x(index2);


% Compute the full width, half max
LSFHPW = x2 - x1;
