function [RadDistH, RadDistV, raw]=radseg(RDseg, namesIndex, path, Labels, Ang)
% RADSEG divides the MTF_Results.mat (the NS-SFRs) into segments according 
% to the position of the ROI in the frame. 
%
% In addtion the LSF Half Peak Width is measured from the resampled ESF
%
% INPUT:
%   RDseg           -       The number of radial segments ('Dohnuts')
%   namesIndex      -       The NS-SFR data to load
%   Path            -       Directory to the data
%   Lables          -       Array that is used to partition the database
%   Ang             -       System e-SFR Estimaion Angle parameter range 
%
% OUTPUT:
%   RadDistH        -       Cell array cotraining the divided NS-SFR data,
%                           edge perameters and the LSF Half Peak Width
%                           for the Horizontal edges
%   RadDistV        -       Cell array cotraining the divided NS-SFR data,
%                           edge perameters and the LSF Half Peak Width
%                           for the Vertical edges
%   raw             -       Is RAW or TIFF? raw==1 or ==0
%
% O. van Zwanenberg (2022)
% 
% UNIVERSITY OF WESTMINSTER 
%              - COMPUTATIONAL VISION AND IMAGING TECHNOLOGY RESEARCH GROUP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Determin if RAW or TIFF
Load_Result=[path '/' namesIndex{1,1}];
load(Load_Result, 'MTF_Results');
if size(MTF_Results,1)==1
    rgbC = 1;
    raw = 0;
else
    rgbC = 3;
    raw = 1;
end

% Load NS-SFR data into cell array - sedmented into the Radial Distances
RadDist=cell(rgbC,RDseg);
data=zeros(1,7);
data=num2cell(data);
for RGB = 1:rgbC
    for a=1:RDseg
        RadDist{RGB,a}=data;
    end
end
RadDistH=RadDist;
RadDistV=RadDist;

ZeroindxH=cell(rgbC,RDseg);
ZeroindxV=ZeroindxH;

% Horizontal and Vertical counter
hv=zeros(rgbC,RDseg);

% RDindex=1/RDseg; 
[Rad, ~] = RadAnnuli(MTF_Results{1, 2} , RDseg);

for A=1:size(namesIndex,1)
    % If the partition label is true, use image
    if Labels(A)==1
        nameIndex=namesIndex(A,1);
        Load_Result=[path '/' nameIndex{1,1}];
%         pause(0.5);
%         MTF_Results = parload(Load_Result, 'MTF_Results');
        load(Load_Result, 'MTF_Results');
        % Orentation
        for O=4:5
            switch O
                case 4 % Horizontal
                    RadDist=RadDistH;
                    Zeroindx = ZeroindxH;
                    HV=hv;
                case 5 % Vertical 
                    RadDist=RadDistV;
                    Zeroindx = ZeroindxV;
                    HV=hv;
            end
            for RGB = 1:size(MTF_Results,1)
                if RGB~=3 && RGB~=4
                    rgb = RGB;
                elseif RGB==3
                    rgb = 2; % Double Green Channel
                elseif RGB==4
                    rgb = 3;
                end
                for B=1:size(MTF_Results{RGB, O},1)
                    if ~isempty(MTF_Results{RGB, O}{B, 1})
                        if MTF_Results{RGB, O}{B, 7}==0
                            if abs(MTF_Results{RGB, O}{B, 6})<=Ang(1,2) && abs(MTF_Results{RGB, O}{B, 6})>=Ang(1,1) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                RD=MTF_Results{RGB, O}{B, 3};
                                a=0;
%                                 b=RDindex;
                                for C=1:RDseg
                                    r = Rad(C);
                                    if RD>=a && RD<r
                                        HV(rgb,C)=HV(rgb,C)+1;
                                        RadDist{rgb, C}{HV(rgb,C),1}=MTF_Results{RGB, O}{B, 1}; % NS-SFR
                                        RadDist{rgb, C}{HV(rgb,C),2}=MTF_Results{RGB, O}{B, 2}; % ESF
                                        RadDist{rgb, C}{HV(rgb,C),3}=MTF_Results{RGB, O}{B, 3}; % Location RD
                                        RadDist{rgb, C}{HV(rgb,C),4}=MTF_Results{RGB, O}{B, 4}; % ESF FWHM
                                        RadDist{rgb, C}{HV(rgb,C),5}=MTF_Results{RGB, O}{B, 5}; % Edge Contrast
                                        RadDist{rgb, C}{HV(rgb,C),6}=MTF_Results{RGB, O}{B, 6}; % Edge Angle
                                        RadDist{rgb, C}{HV(rgb,C),7}=[];                        % empty
                                        % Check for zero ESF widths & remove - due to NAN grey 
                                        % SFR or flat grey ESF
                                        if MTF_Results{RGB, O}{B, 4}==0
                                            Zeroindx{rgb,C}{HV(rgb,C),1}=1;
                                        else 
                                            Zeroindx{rgb,C}{HV(rgb,C),1}=0;
                                        end
                                a=r;
%                                     b=b+RDindex;
                                    end 
                                end 
                            end
                        end  
                    end
                end
            end
            switch O
                case 4 % Horizontal
                    RadDistH= RadDist;
                    ZeroindxH = Zeroindx;
                    HV=hv;
                case 5 % Vertical 
                    RadDistV= RadDist;
                    ZeroindxV = Zeroindx;
            end
        end
    end
end
% Remove erros
for RGB = 1:rgbC
    for A=1:RDseg
        if isempty(ZeroindxH{RGB,A})==0
            err=find([ZeroindxH{RGB,A}{:}]==1);
            emp=isempty(err);
            if emp==0
                ZeroindxH{RGB,A}(err,:)=[];
                RadDistH{RGB,A}(err,:)=[];
            end
        end
        if isempty(ZeroindxV{RGB,A})==0
            err=find([ZeroindxV{RGB,A}{:}]==1);
            emp=isempty(err);
            if emp==0
                ZeroindxV{RGB,A}(err,:)=[];
                RadDistV{RGB,A}(err,:)=[];
            end
        end
    end
end

