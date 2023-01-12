function [BWh, BWv, Mh, Mv, SEh, SEv]=findEdge(img, ST, esfW, MinEdge)
% findEdge is an adaption of the canny edge detector, providing vertical
% and horizontal edges sepratly, edge masks and determines step-edges.
%
% Input:
%    img:       The Image.
%    ST:        (Optional) Adjusible varible is the isStepEdege 
%               threshold, this  value should be 0-1 (defult is 0.02 
%               normialised pixel  value, aproximaly equal to +/-5 pixel 
%               value for 8 bit image).
%   esfW:       (Optional) The ESF width, based on the MTF model (defult is
%               5 pixels). 
%   MinEdge:    (Optional) Minimum edge height - defult 20 pixels
%
% Output:
%    BWh:       Horizontal Edge Location Binary Image
%    BWv:       Vertical Edge Location Binary Image
%    Mh:        Horizontal Edge Magnitude Image
%    Mv:        Vertical Edge Magnitude Image
%
% Adpated from MATLAB function [BW,thresh] = EDGE(I,'canny',...)
% Adapted by O. van Zwanenberg 2022, University of Westminster

%--------------------------------------------------------------------------
% Convert img to grayscale and double
if size (img, 3)>1
    img = rgb2gray(img);
end
% convert to normalised double
img=im2double(img);

% Value for Thresholding
% if ~exist('Con','var')
%     Con = [0,1];
% end
if ~exist('ST','var')
    ST = 0.04;
end
if ~exist('esfW','var')
    esfW = 5;
end
if ~exist('MinEdge','var')
    MinEdge = 20;
end

% Magic numbers
PercentOfPixelsNotEdges = .7; % Used for selecting thresholds
ThresholdRatio = .4;          % Low thresh is this fraction of the high.
% Calculate gradients using a derivative of Gaussian filter  
% Calculate directions/orientations
% [dx, dy] = smoothGradient(img, sqrt(2)); %defult sigma = sqrt(2)
% Gaussian Filter Coefficient
B = [2, 4, 5, 4, 2; 4, 9, 12, 9, 4;5, 12, 15, 12, 5; ...
    4, 9, 12, 9, 4;2, 4, 5, 4, 2 ];
B = 1/159.* B;

% Convolution of image by Gaussian Coefficient
A=conv2(img, B, 'same');

% Filter for horizontal and vertical direction
KGv = [-1, 0, 1; -2, 0, 2; -1, 0, 1];
KGh = [1, 2, 1; 0, 0, 0; -1, -2, -1];
% Convolution by image by horizontal and vertical filter
dx = conv2(A, KGv, 'same');
dy = conv2(A, KGh, 'same');


arah = atan2 (dy, dx);
arah = arah*180/pi;

pan=size(A,1);
leb=size(A,2);

% Adjustment for negative directions, making all directions positive
IJ = find(arah<0);
arah(IJ)=360+arah(IJ);

DX=zeros(pan, leb);
DY=zeros(pan, leb);

% Split Horizontal and Vertical according to measure edge direction
IJ = find(arah>= 0  & arah < 45 | ...
                    arah>= 135  & arah < 225 |...
                    arah>= 315  & arah <= 360);
DX(IJ) = dx(IJ);
IJ = find(arah>= 45  & arah < 135 |...
                    arah>= 225  & arah < 315);
DY(IJ) = dy(IJ);

% Orentation -  H and then V
for O=1:2
    switch O
        case 1
            % Calculate Magnitude of Gradient
            magGrad = hypot(zeros(size(DY)), DY);
        case 2
            magGrad = hypot(DX, zeros(size(DX)));
    end
    % Normalize for threshold selection
    magmax = max(magGrad(:));
    if magmax > 0
        magGrad = magGrad / magmax;
    end

    % Determine Hysteresis Thresholds
    [lowThresh, highThresh] = selectThresholds([], magGrad, ...
        PercentOfPixelsNotEdges, ThresholdRatio, mfilename);

    % Perform Non-Maximum Suppression Thining and Hysteresis Thresholding 
    % of Edge Strength
    e = thinAndThreshold(DX, DY, magGrad, lowThresh, highThresh);
    switch O
        case 1 
            e = rot90(e);
            magGrad = rot90(magGrad);
            bwh = e;
%             Mh = magGrad;
        case 2
            bwv = e;
%             Mv = magGrad;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Cull Unwanted Edges %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remove all edges smaller then 20 pixels in length:
bwh= EdgeReduction(bwh, MinEdge);
bwv= EdgeReduction(bwv, MinEdge);

% Convert the binary arrays to list of cordinates [x,y]
[BWcoh]=BW2Co(bwh);
[BWcov]=BW2Co(bwv);

% Proximity Filter, removes detected edges that are too close
% i.e., remove textures - Set using esfW
bwh=Proximity(bwh,BWcoh, 1, esfW);
bwv=Proximity(bwv,BWcov, 1, esfW);
% Re-Apply: Remove all edges smaller then 20 pixels in length:
bwh= EdgeReduction(bwh, MinEdge);
bwv= EdgeReduction(bwv, MinEdge);
%--------------------------------------------------------------------------
% Contrast thresholds, step-edge determination and edge of intrest
% Masking
[ImGradV, ImGradH]=gradient(img);
ImGradH = rot90(ImGradH);

% Thrshold the Gradient Array - anything between -Var2 and Var2 = 0,
% else if positve = 1 or negative =-1
SE = strel('square',2);
SE2 = strel('square',7);
for O=1:2
    switch O
        case 1 % H
            ImGrad = ImGradH; %./max(max(ImGradH));
            bw = bwh;
%             im = rot90(img);
            stepedge = zeros(size(bw));
        case 2 % V
            ImGrad = ImGradV; %./max(max(ImGradV));
            bw = bwv;
%             im = img;
            stepedge = zeros(size(bw));
    end

    BW = zeros(size(bw));

    ImGrad(ImGrad<ST & ImGrad>-ST) = 0;
    ImGrad(ImGrad<0) = -1;
    ImGrad(ImGrad>0) = 1;

    [co]=BW2Co(bw);
    bwmask = bwselect(ImGrad,co(:,1),co(:,2),4);
    Mask = zeros(size(bwmask));

    for a = 1:size(bw, 1)
        Mask(a,:) = imreconstruct(bwmask(a,:) & bw(a,:), ...
            bwmask(a,:) | bw(a,:), 4);
    end
    Mask = imerode(Mask, SE);
    Mask = imdilate(Mask, SE2);
    Mask =logical(Mask);
    MaskGrad = zeros(size(Mask));
    MaskGrad(Mask)= ImGrad(Mask);
%     Mask2=zeros(size(Mask));
%     BW2 = Mask2;
    [Y,~,iy] = unique(co(:,2));

    for yy = 1:length(Y)
        X=co(iy==yy,1);
%         y=zeros(size(X));
%         y(:)=Y(yy);
        labeledMask = bwlabel(Mask(Y(yy),:));
        for ii = 1:max(labeledMask)
            % Is the edge location central to the Mask? 
            Indx = find(labeledMask==ii);
            xx=find(ismember(X, Indx)==1);
            if size(X(xx),1)>1 || isempty(xx)
                continue
            end
            if X(xx)== Indx(1,1) || X(xx)== Indx(1,end)
                 continue
            end
            
            % Is Step Edge? (measures the gradent of a profile to determine
            % if there is a single positive or negarive gradient,
            % otherwise stepedge=0)
            gpos = 0;
            gneg = 0;
            grad = MaskGrad(Y(yy),Indx);
            [N, ~] = find(grad==-1);
            [P, ~] = find(grad==1);
            N=sum(N);
            P=sum(P);

            if N<P
                gpos = 1;
                % Only egdes of above gradient threshold in BW
                BW(Y(yy), X(xx)) = 1;
            elseif P<N
                gneg = 1;
                % Only egdes of above gradient threshold in BW
                BW(Y(yy), X(xx)) = 1;
            end
            % Create SE Map
            if gpos==1 && gneg==0
                stepedge(Y(yy),X(xx))=1;
            elseif gneg==1 && gpos==0
                stepedge(Y(yy),X(xx))=-1;
            else
                stepedge(Y(yy),X(xx))=0;
            end
            
% % %             % If step edge is detected measure contrast
% % %             if stepedge == 1
% % %                 % Do not apply contrast limits if [0,1], measured with indx
% % %                 % with padding
% % %                 if  Indx(1,1) == 1
% % %                     Ipad1 = Indx(1,1);
% % %                 else
% % %                     Ipad1 = Indx(1,1)-1;
% % %                 end
% % %                 if Indx(1,end) == size(im,2)
% % %                     Ipad2 = Indx(1,end);
% % %                 else
% % %                     Ipad2 = Indx(1,end)+1;
% % %                 end
% % %                 if Con(1,1)~=0 || Con(1,2)~=1
% % %                     i1 = mean(im(Y(yy), Ipad1:X(xx)-1));
% % %                     i2 = mean(im(Y(yy), X(xx)+1:Ipad2));
% % %                     con = abs(i1-i2)/(i1+i2);
% % %                     con = floor(con) + ceil( (con-floor(con))/0.05) * 0.05;
% % %                     if con<=Con(1,2) && con>=Con(1,1)
% % %                         C = 1;
% % %                     else
% % %                         C = 0;
% % %                     end
% % %                 else
% % %                     C = 1;
% % %                 end
% % %             
% % %                 if C==1
% % %                     BW2(Y(yy), X(xx)) = 1;
% % %                     if gpos==1
% % %                     Mask2(Y(yy), Ipad1:Ipad2) = 1;
% % %                     elseif gneg==1
% % %                     Mask2(Y(yy), Ipad1:Ipad2) = -1;
% % %                     end
% % %                 end
% % %             end
%             MaskArea = bwselect(Mask(Y(yy),:),X(xx), 1,4);
%         Row=MaskGrad(Y(yy),:);
        end
    end
    switch O
    case 1 % H
        BWh = logical(BW);
        Mh =  Mask; 
        SEh = stepedge;
    case 2 % V
        BWv = logical(BW);
        Mv =  Mask; 
        SEv = stepedge;
    end
end

% % % % Split positve and negative Step edge prfiles
% % % BWh = cell(1,2);
% % % BWv = BWh;
% % % bh = zeros(size(bwh));
% % % bv = zeros(size(bwv));
% % % BWh{1,1}=bh;
% % % BWh{1,2}=bh;
% % % BWv{1,1}=bv;
% % % BWv{1,2}=bv;
% % % 
% % % Mh  = BWh;
% % % Mv  = BWv;
% % % 
% % % % Mh{1,1}= imreconstruct(mh==1 & bwh, mh==1 | bwh, 4);
% % % % Mh{1,2}= imreconstruct(mh==-1 & bwh, mh==-1 | bwh, 4);
% % % % Mv{1,1}= imreconstruct(mv==1 & bwv, mv==1 | bwv, 4);
% % % % Mv{1,2}= imreconstruct(mv==-1 & bwv, mv==-1 | bwv, 4);
% % % for Y = 1:size(bwh,1)
% % %     Mh{1,1}(Y,:) = imreconstruct(mh(Y,:)==1 & bwh(Y,:), ...
% % %         mh(Y,:)==1 | bwh(Y,:), 4);
% % %     Mh{1,2}(Y,:) = imreconstruct(mh(Y,:)==-1 & bwh(Y,:), ...
% % %         mh(Y,:)==-1 | bwh(Y,:), 4);
% % % end
% % % for Y = 1:size(bwv,1)
% % %     Mv{1,1}(Y,:) = imreconstruct(mv(Y,:)==1 & bwv(Y,:), ...
% % %         mv(Y,:)==1 | bwv(Y,:), 4);
% % %     Mv{1,2}(Y,:) = imreconstruct(mv(Y,:)==-1 & bwv(Y,:), ...
% % %         mv(Y,:)==-1 | bwv(Y,:), 4);
% % % end
% % % BWh{1,1}= bwh & Mh{1,1};
% % % BWh{1,2}= bwh & Mh{1,2};
% % % BWv{1,1}= bwv & Mv{1,1};
% % % BWv{1,2}= bwv & Mv{1,2};
% % % % Remove all edges smaller then 20 pixels in length:
% % % BWh{1,1} = EdgeReduction(BWh{1,1}, 20);
% % % BWh{1,2} = EdgeReduction(BWh{1,2}, 20);
% % % BWv{1,1} = EdgeReduction(BWv{1,1}, 20);
% % % BWv{1,2} = EdgeReduction(BWv{1,2}, 20);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Local Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %   smoothGradient
% function [GX, GY] = smoothGradient(I, sigma)
% 
% % Create an even-length 1-D separable Derivative of Gaussian filter
% 
% % Determine filter length
% filterExtent = ceil(4*sigma);
% x = -filterExtent:filterExtent;
% 
% % Create 1-D Gaussian Kernel
% c = 1/(sqrt(2*pi)*sigma);
% gaussKernel = c * exp(-(x.^2)/(2*sigma^2));
% 
% % Normalize to ensure kernel sums to one
% gaussKernel = gaussKernel/sum(gaussKernel);
% 
% % Create 1-D Derivative of Gaussian Kernel
% derivGaussKernel = gradient(gaussKernel);
% 
% % Normalize to ensure kernel sums to zero
% negVals = derivGaussKernel < 0;
% posVals = derivGaussKernel > 0;
% derivGaussKernel(posVals) = derivGaussKernel(posVals)/sum(derivGaussKernel(posVals));
% derivGaussKernel(negVals) = derivGaussKernel(negVals)/abs(sum(derivGaussKernel(negVals)));
% 
% % Compute smoothed numerical gradient of image I along x (horizontal)
% % direction. GX corresponds to dG/dx, where G is the Gaussian Smoothed
% % version of image I.
% GX = imfilter(I, gaussKernel', 'conv', 'replicate');
% GX = imfilter(GX, derivGaussKernel, 'conv', 'replicate');
% 
% % Compute smoothed numerical gradient of image I along y (vertical)
% % direction. GY corresponds to dG/dy, where G is the Gaussian Smoothed
% % version of image I.
% GY = imfilter(I, gaussKernel, 'conv', 'replicate');
% GY  = imfilter(GY, derivGaussKernel', 'conv', 'replicate');
%--------------------------------------------------------------------------
%  selectThresholds
function [lowThresh, highThresh] = selectThresholds(thresh, magGrad, ...
    PercentOfPixelsNotEdges, ThresholdRatio, ~)

[m,n] = size(magGrad);

% Select the thresholds
if isempty(thresh)
    counts=imhist(magGrad, 64);
    highThresh = find(cumsum(counts) > PercentOfPixelsNotEdges*m*n,...
        1,'first') / 64;
    lowThresh = ThresholdRatio*highThresh;
elseif length(thresh)==1
    highThresh = thresh;
    if thresh>=1
        error(message('images:edge:singleThresholdOutOfRange'))
    end
    lowThresh = ThresholdRatio*thresh;
elseif length(thresh)==2
    lowThresh = thresh(1);
    highThresh = thresh(2);
    if (lowThresh >= highThresh) || (highThresh >= 1)
        error(message('images:edge:thresholdOutOfRange'))
    end
end
%--------------------------------------------------------------------------
%   thinAndThreshold
function H = thinAndThreshold(dx, dy, magGrad, lowThresh, highThresh)
% Perform Non-Maximum Suppression Thining and Hysteresis Thresholding of
% Edge Strength

% We will accrue indices which specify ON pixels in strong edgemap
% The array e will become the weak edge map.

E = images.internal.builtins.cannyFindLocalMaxima(dx,dy,magGrad,lowThresh);

if ~isempty(E)
    [rstrong,cstrong] = find(magGrad>highThresh & E);

    if ~isempty(rstrong) % result is all zeros if idxStrong is empty
        H = bwselect(E, cstrong, rstrong, 8);
    else
        H = false(size(E));
    end
else
    H = false(size(E));
end