% E-SFR Estimation from NS-SFR Data
% 
% Copyright (c) 2023 O. van Zwanenberg
% UNIVERSITY OF W1ESTMINSTER PhD Reserch
%              - COMPUTATIONAL VISION AND IMAGING TECHNOLOGY RESEARCH GROUP
% Director of Studies:  S. Triantaphillidou
% Supervisory Team:     R. Jenkin & A. Psarrou

clc; close all; 

Ang=[2.5, 42.5]; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PERAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Number of Radial Distance Segemnts ('Radial Distance Dohnuts')
RDseg=3;
% Percentile of the LSF half peak width distribution to be used in the 
% sys-SFR estimation
Percentile = 20; % Top 10th percentile of sharpest edges per RDseg
Percentile2=[];

% Thresholds
% System e-SFR Estimaion parameter ranges
sysPar = [2,      35;     0.55,   0.65;   20,      130    ];
%        [minAng, maxAng; minCon, maxCon; minROIh, maxROIh];
% SFR mean Weights
AveW = [1.00  , 0.75    , 0.5    ];
%      [Centre, Part-Way, Corners];  - [1] = no weight
%                                    - [1,2,3,...,n] = more radial regions  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% OPEN DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mat = uigetdir([], 'Select folder conatining NS-SFR Data');
matemp=isempty(mat);
if matemp==0
    load([mat '/ImageNamesIndex.mat']); 
    nameIndex=namesIndex(1,1);
    Load_Result=[mat '/' nameIndex{1,1}];
    load(Load_Result); 
else
    disp('No Data in chosen Folder');
    return 
end
tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% SEGMENT NS-SFRS BY RD %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% No Label - Use to partition the NS-SFR data acording to scene types
label = zeros(size(namesIndex,1),1)+1;
% Segment
[RadDistH, RadDistV, raw] =radseg(RDseg, namesIndex, mat, label, Ang);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LSF FWHM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Obtain the strongest NS-SFRs    
% Take a distribution of the range of LSF FWHM
% Measure Area and hight of the NS-SFR lobes
for RGB=1:size(RadDistH,1)
    % Horizontal Edges
    RadDistH(RGB,:)=nssfrScore(RadDistH(RGB,:),RDseg, Percentile,...
        Percentile2);
    % Vertical Edges
    RadDistV(RGB,:)=nssfrScore(RadDistV(RGB,:), RDseg, Percentile, ...
        Percentile2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% MEADIAN NS-SFR PER RD %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LSFH = cell(1, RDseg);
LSFV = LSFH;
for RGB=1:size(RadDistH,1)
    for A=1:RDseg
        a = cell2mat(RadDistH{RGB, A}(:,7));
        if sum(a)==0
            continue
        end
        s = find(a==1)';
        LSFH{RGB,A}=RadDistH{RGB, A}((all((a==1),2)),1:2);
        
        a = cell2mat(RadDistV{RGB, A}(:,7));
        if sum(a)==0
            continue
        end
        LSFV{RGB,A}=RadDistV{RGB, A}((all((a==1),2)),1:2);
    end
end

[AveSFRh, AveLSFh]=medianSFR(LSFH, AveW, raw);

[AveSFRv, AveLSFv]=medianSFR(LSFV, AveW, raw);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SAVE DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save([mat '/eSFR_Estimation.mat'], 'AveSFRh', 'AveLSFh', ...
    'AveSFRv', 'AveLSFv', 'RadDistH', 'RadDistV')
toc
disp('Compleated System e-SFR Estimation');
% clearvars 
beep
