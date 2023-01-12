function [AveSFR, AveLSF]=medianSFR(LSF, RD, RAW)
% medianSFR takes a mean Spatial Frequency Response (SFR) through  
% registering the LSFs (from the supersampled ESFs) and then taking the 
% modulas Fourier Transform of the median LSF
%
% The SFR has been caculated using code from sfrmat4  
% (Copyright (c) Burns Digital Imaging, 2020. [Online]. 
% Available: http://burnsdigitalimaging.com/software/sfrmat/.)
% 
%
% Input: 
%       LSF       -    Cell Array of Edge Spread Functions (LSFs) to be 
%                      used and their original measured NS-SFR. 
%                      The array should have column 2 being LSFs and column
%                      1 being the radial distance. 
%       RD        -    This the radial distance array. The number of
%                      elements in the array refers to the number of radial
%                      distance segments. Each element should be the
%                      weighting of the mean for that radial distance
%                      segment, e.g. splitting the frame into three parts,
%                      centure, partway and cornors with weightings of
%                      1.00, 0.75 and 0.50 respectively would be shown as
%                      [1, 0.75, 0.5]. The defult is RD=[1], i.e. no radial
%                      distance segmentation and no weighting
%       RAW       -    Is the data RAW, True=1, False=0. If True frequency
%                      is halfed and a SFR for each channel is calculated. 
% Output: 
%       avesfr    -    The output mean SFR
%       AveLSF    -    The output mean LSF
%
% 2022, O. van Zwanenberg
% UNIVERSITY OF WESTMINSTER 
%              - COMPUTATIONAL VISION AND IMAGING TECHNOLOGY RESEARCH GROUP

%--------------------------------------------------------------------------
switch nargin
    case 1
        RD=1; % Defult RD
        RAW = 0;
    case 2
        if size(RD,2)==1 && RD~=1
          disp(['If there is only one Radial Segment (RD), ' ...
              'then the weighting value should equal 1.00, i.e. RD=[1]']);
        return   
        end
        RAW =0;
    case 3
        if size(RD,2)==1 && RD~=1
          disp('If there is only one Radial Segment (RD) then the weighting value should equal 1.00, i.e. RD=[1]');
        return   
        end
    otherwise
        disp('Incorrect number or arguments');
        return 
end
%--------------------------------------------------------------------------
if RAW==1
    RGB=3;
else
    RGB=1;
end

AveSFR=cell(2,size(LSF, 2)+1);
AveLSF=cell(2,size(LSF, 2)+1);

for rgb=1:RGB
    % Find Mid LSF
    L = 0;
    % Index=zeros(1,1);
    % a=0; 
    Mid = cell(size(LSF,2),1);
    maxMid = zeros (size(LSF,2), 1); 
    for A =1:size(LSF,2)
        if isempty (LSF{rgb,A})
            Mid{A,1}=0;
            maxMid(A,1) = 0;
            continue
        end
        mid = zeros (size(LSF{rgb, A},1), 1); 
        for B=1:size(LSF{rgb, A},1)
            mi = find(abs(LSF{rgb,A}{B,2})==max(abs(LSF{rgb,A}{B,2})));
            mid(B,1)=mi(1,1); % first mid
            l=length(LSF{rgb, A}{B,2});
            if l>L
                L=l;
            end
        end
        Mid{A,1}=mid;
        maxMid(A,1) = max(mid);
    end
    maxMid=max(maxMid);
    % Register the LSFs
    LSFrd=cell(1,size(LSF,2)+1);
    
    for A=1:size(LSF,2)
        if Mid{A,1}==0
            LSFrd{1,A}=[];
            continue
        end
        LSFreg=zeros(size(LSF{rgb,A},1),L+maxMid);
        LSFrd{1,A}=LSFreg;
    end
    % Frame Radial distances
    for A=1:size(LSF,2)  
        if Mid{A,1}==0
            continue
        end
    %     l=0;
        for B=1:size(LSF{rgb,A},1)
    %         l=l+1;
            m=Mid{A,1}(B,1);
            dif=maxMid-m;
            if dif~=0
                ld=length(LSF{rgb,A}{B,2})+(dif-1);
            else
                ld=length(LSF{rgb,A}{B,2});
                dif=1;
            end
            % Normalise the Area of the LSFs
            lsf = LSF{rgb,A}{B,2};
            Area=trapz(1:size(lsf,2),lsf);
            normLSF=lsf./Area;
            LSFrd{1,A}(B,dif:ld)=normLSF;
        end
    end
    
    % Measure the median LSF per Radial Distance
    LSFRD=zeros(size(LSF,2),L+maxMid);
    LSFw = zeros(size(LSFRD));
    RDw=zeros(size(LSF,2),1);
    i1=1;
    i=size(LSF,2)/size(RD,2);
    i2=i;
    for A=1:size(RD,2)
        RDw(i1:i2,1)=RD(1,A);
        i1=i1+i;
        i2=i2+i;
    end
    
%     AveLSF = cell(size(LSFrd));
    for A=1:size(LSF,2) 
        if size (LSFrd{1,A},1)>1
            LSFRD(A,:) = median(LSFrd{1,A});
        elseif size (LSFrd{1,A},1)==1
            LSFRD(A,:) = LSFrd{1,A};
        else
            LSFRD(A,:) = NaN;
        end
        AveLSF{rgb,A}=LSFRD(A,:);
        LSFw(A,:)  = RDw(A);
    end
    
    % Remove NaNs
%     nans=isnan(LSFRD);
    LSFw(sum(isnan(LSFRD), 2) ~= 0, :) = [];
    LSFRD(sum(isnan(LSFRD), 2) ~= 0, :) = [];
    
    % Take the weighted median LSF across all ditance segments
    AveLSF{rgb,size(LSF, 2)+1} = zeros(1,L+maxMid);
    for A=1:L+maxMid
         AveLSF{rgb, size(LSF, 2)+1}(A) = weightedMedian(LSFRD(:,A),LSFw(:,A));
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     AveSFR = cell(1, size(LSF, 2)+1);
    % Take the SFR from the mean LSF (AveLSF) 
    %                - Adapted from sfrmat4 (Copyright (c) 2020 Peter 
    %                  D. Burns)
    
    % % % delfac = cos(atan(vslope));    
    % % % del = 1*delfac;  % input pixel sampling normal to the edge
    % % % del2 = del/4;   % Supersampling interval normal to edge
    del2=1/4; % (del/nbin)
    nn= L+maxMid;% floor(npix *nbin);
    
    nn2 =  floor(nn/2) + 1;
    % frequency 
    freq = zeros(nn, 1);
    for n=1:nn   
        freq(n) = (n-1)/(del2*nn);
    end
    if RAW==1
        freq=freq/2;
    end
    % [correct] = fir2fix(n, m);
    % Correction for MTF of derivative (difference) filter
     % dcorr corrects SFR for response of FIR filter
    dcorr = ones(nn2, 1);
    m=3-1;
    scale = 1;
    for i = 2:nn2
        dcorr(i) = abs((pi*i*m/(2*(nn2+1))) / sin(pi*i*m/(2*(nn2+1))));
        dcorr(i) = 1 + scale*(dcorr(i)-1);
      if dcorr(i) > 10  % Note limiting the correction to the range [1, 10]
        dcorr(i) = 10;
      end
    end
    % desired frequency intervals
    uq=(0:0.01:0.50)';
    for A=1:size(LSF, 2)+1
        % if only one edge, use original NS-SFR measument
        if size (LSFrd{1,A},1)==1
            AveSFR{rgb, A} = LSF{1,A}{1,1};
            continue
        end
        if ~isempty (AveLSF{1,A}) && ~isnan (AveLSF{1,A}(1,1))
            sfr =  zeros(nn, 1);
            temp = abs(fft(AveLSF{1,A}, nn));
            sfr(1:nn2, 1) = temp(1:nn2)/temp(1);
            sfr(1:nn2, 1) = sfr(1:nn2, 1).*dcorr(1:nn2);
            if ~isempty(sfr)
                mq=interp1(freq, sfr, uq, 'pchip');
                AveSFR{rgb, A}(:,1)=uq;
                AveSFR{rgb, A}(:,2)=mq;
            else
                AveSFR{rgb, A}=[];
            end
        else 
            AveSFR{rgb, A}=[];
        end  
    end
end