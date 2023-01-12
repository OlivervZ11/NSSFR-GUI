function [R, RD] = RadAnnuli(im , num)
% RadAnnuli - Segments the image frame into even area Radial Annuli
% Input:
%   im       Image frame to segment
%   num      The number if Radial Anuli to segment the frame into
%
% Output:
%   R        The Radius of each Anuli
%   RD       Segmented frame 
%   
% Copyright (c) 2021 O. van Zwanenberg
% UNIVERSITY OF W1ESTMINSTER PhD Reserch
%              - COMPUTATIONAL VISION AND IMAGING TECHNOLOGY RESEARCH GROUP

% Image area per segment
[Y, X, ~] = size(im);
% im centure
y = round(Y/2);
x = round(X/2); 

%Measure the disance from centure to corner, i.e. max distance 
CO = [y,x;Y,X];
MAXdist = pdist(CO,'euclidean');

% A = pi*MAXdist^2;
A =  pi*x^2;
a = A/(num);

% Set up Radus list
R = zeros(1,num);

% Loop num to caculate R
for N = 1:num
    r = sqrt((a*N)/pi);
    R(1,N)=r;
% R(1,N)=(MAXdist/num)*N;
%     if r>Y
%         % Caulate Area of Chord
%         ang = 2*acos(a-((r-Y)/r));
%         a2 = 0.5*(R^2)*(ang-sin(ang));
%         A2=a2*2;
%         % Add area of chordx2  
%         % - Will not work, the area is distrubuted throughout the circle
%         a3 = a+A2;
%         R(1,N) = sqrt(a3/pi);
%     else
%         R(1,N)=r;
%     end
end
%--------------------------------------------------------------------------
RD = zeros(Y, X); 
[columnsInImage, rowsInImage] = meshgrid(1:X, 1:Y);
R(1,num)=MAXdist;
for rd = num:-1:1
    radius = R(1,rd);
    circlePixels = (rowsInImage - y).^2 ...
    + (columnsInImage - x).^2 <= radius.^2;
    RD(circlePixels)=rd;
end
% % XY = round(MAXdist*2);
% % xy = round(MAXdist);
% % RD = zeros(XY, XY); 
% % [columnsInImage, rowsInImage] = meshgrid(1:(XY), 1:(XY));
% % for rd = num:-1:1
% %     radius = R(1,rd);
% %     circlePixels = (rowsInImage - xy).^2 ...
% %     + (columnsInImage - xy).^2 <= radius.^2;
% %     RD(circlePixels)=rd;
% % end
 I=zeros(1,num);
for rd = 1:num
    [~,i] = find (RD ==rd);
    I(1,rd)=size(i,1);
end
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Correct the loss of area 

i = I-(mean(I));
aa=0;
for N = 1:num
    r = sqrt(((a-i(N))+aa)/pi);
    R(1,N)=r;
    aa=aa+(a-i(N));
end
%--------------------------------------------------------------------------
RD = zeros(Y, X); 
[columnsInImage, rowsInImage] = meshgrid(1:X, 1:Y);
R(1,num)=MAXdist;
for rd = num:-1:1
    radius = R(1,rd);
    circlePixels = (rowsInImage - y).^2 ...
    + (columnsInImage - x).^2 <= radius.^2;
    RD(circlePixels)=rd;
end

 I=zeros(1,num);
for rd = 1:num
    [~,i] = find (RD ==rd);
%     disp(['Area of Radial Annuli ' num2str(rd) ' = ' num2str(size(i,1))]);
    I(1,rd)=size(i,1);
end


R=R/MAXdist;