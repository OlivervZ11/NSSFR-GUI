% Natural Scenes derived Spatial Frequency Response (NS-SFR) Extraction
% 
% Copyright (c) 2023 O. van Zwanenberg
% UNIVERSITY OF W1ESTMINSTER
%              - COMPUTATIONAL VISION AND IMAGING TECHNOLOGY RESEARCH GROUP
% Director of Studies:  S. Triantaphillidou
% Supervisory Team:     R. Jenkin & A. Psarrou

clc; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Varibles
% Contrast Range 
Con=[0.55, 0.65]; 
% isStepEdge zeros gradient value 
ST   = 0.02;
% ESF min width for system (pixels) - Change this threshold based on
% a modeled MTF - pixels for no overlapping ESFs
esfW = 5; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% ----- READ FILES ----- %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RAW(.dng) or TIFF(.tif) or PNG(.png) or JPEG(.jpg)? 
answer = questdlg('Image Dataset File Format', ...
	'Data', ...
	'RAW(.dng)','Other', 'RAW(.dng)');
switch answer
    case 'RAW(.dng)' 
        raw = 1;
    case 'Other'
        raw = 0; 
end
if raw==0
    answer = questdlg('Image Dataset File Format', ...
	    'Data', ...
	    'TIFF(.tif)', 'PNG(.png)', 'JPEG(.jpg)', 'TIFF(.tif)');
    switch answer
        case 'TIFF(.tif)'
            ft = 1;
        case 'PNG(.png)'
            ft = 2;
        case 'JPEG(.jpg)'
            ft = 3;
    end
end

% Read image file names from user Folder
selpath = uigetdir([], 'Select folder conatining image dataset');
if raw == 1
    imfiles=dir(fullfile([selpath '\*.dng']));
else
    if ft ==1
        imfiles=dir(fullfile([selpath '\*.tif']));
    elseif ft == 2
        imfiles=dir(fullfile([selpath '\*.png']));
    elseif ft == 3
        imfiles=dir(fullfile([selpath '\*.jpg']));
    end
end

% imnumber stores the number of files that have been read
imnumber=size(imfiles,1);
if imnumber==0
   disp('No Images found in selected folder');
   beep
   return
end
disp(['Number of Detected images = ' num2str(imnumber)]);

% Select folder to save NS-SFR data
resultdir = uigetdir(selpath, 'Select folder to save NS-SFR data');

% If NS-SFR Data already exisits in dir, code will continue to add to the
% file
mat = dir([resultdir '/*.mat']); 
matemp=isempty(mat);
if matemp==0
    load([resultdir '/ImageNamesIndex.mat']); 
    ContinueA=length(namesIndex);
else 
    ContinueA=0;
end

% Large loop of functions to extract edges and MTFs from each image in imfiles structure
time = zeros(imnumber,1);
for A=1:imnumber
    tic
    % Display Progress  
    Waitbartex=['Processing Image...' num2str(A) '/' num2str(imnumber)];
    disp(Waitbartex);
    
    % Read Image A
    if raw == 1
        % Read RAW file - imreadDNG is based on a Nikon D800 .NEF RAW file
        [Ir, Ig1, Ig2, Ib, ~]=imreadDNG(fullfile(imfiles(A).folder,...
            imfiles(A).name), 0);
        % Number of colour chanels to process
        CC = 4; % RGGB
    else
        Im=imread(fullfile(imfiles(A).folder, imfiles(A).name));
        Im=Im(:,:,1:3);    
        % Convert to greyscale 
        im=rgb2gray(Im);
        CC = 1; % Greyscale
    end
    % set up array to save NS-SFRs
    MTF_Results=cell(CC, 5);

    % Loop colour channels (for RAW)
    for B =1:CC 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%% ----- PREP IMAGE ----- %%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % RAW data colour channel switch
        if raw == 1
            switch B
                case 1 % Red Channel
                    Im = Ir;
                    Im_lin=Ir;
                case 2 % Green 1 Channel 
                    Im = Ig1;
                    Im_lin=Ig1;
                case 3 % Green 2 Channel 
                    Im = Ig2;
                    Im_lin=Ig2;
                case 4 % Blue Channel 
                    Im = Ib;
                    Im_lin=Ib;
            end
        end
        %Determine the orientation of the image (if portrait, flip 90deg)
        [x ,y, ~]=size(Im);
        if y<x
           Im=rot90(Im); 
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%% ----- LINEARIZE IMAGE ----- %%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % TIFF files require linerisation
        if raw ==0
            % Extract Image Metadata and Stats
            % ImageInfo=ImInfo(imfiles, A, Im); % NIKON FILE SPESIFCIC -
            % ERROR OCCURS WHEN NO METADATA
            ImageInfo=[];

           % Determine the Colour Space
            if isempty(ImageInfo)==0
                ColourSpace=ImageInfo.ColourSpace;
            else 
                % assume RGB
                ColourSpace = 'RGB';
            end
            % Linearize the ime acording to the Colour Space
            if isequal(ColourSpace ,'sRGB')
                Im_lin=rgb2lin(Im,'ColorSpace', 'sRGB');
            elseif isequal(ColourSpace ,'RGB')
                Im_lin=rgb2lin(Im,'ColorSpace', 'adobe-rgb-1998');
            else
                % assume RGB
                Im_lin=rgb2lin(Im,'ColorSpace', 'adobe-rgb-1998');
            end
        else
             ImageInfo=[];
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%% ----- EDGE DETECTION ----- %%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Apply Edge Detection (based on Canny Edge Detection)
        [BWh, BWv, Mh, Mv, SEh, SEv]=findEdge(Im_lin, ST, esfW);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%% ----- ISOLATE ROIS ----- %%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Isolate, extract and store the ROIs
        ROIsH=ROIisolate(rot90(Im_lin), BWh, Mh, SEh, Con);
        ROIsV=ROIisolate(Im_lin, BWv, Mv, SEv, Con);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%% ----- MEASURE NS-SFR ----- %%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % frequency range
        uq=(0:0.01:0.50)';
        
        % 
        for O=1:2
            switch O
                case 1 
                    R = cell(size(ROIsH,1), 8);
                    r = ROIsH;
%                     rps = PSroiH; 
                    i=rot90(Im_lin);
                case 2
                    R = cell(size(ROIsV,1), 8);
                    r = ROIsV;
%                     rps = PSroiV; 
                    i=Im_lin;
            end
            for a = 1:size(R,1)
                if ~isempty(r{a,1})
                    if size(r{a,1},3)==3
                        ROI=rgb2gray(r{a,1});
                    else
                        ROI = r{a,1};
                    end
%                     ROIps = rps{a,1};
                    mask = r{a,3};
                    %check Canny is at end of ROI
                    yend=find(r{a,2}(end,:)==1, 1);
                    if isempty(yend)
                        ROI=ROI(1:end-1,:);
                        mask=mask(1:end-1,:);
                    end
%                     [~, dat, ~, ~, R{a,2}, ~, ~, ~, R{a,6}, R{a,7}] = ...
%                         sfrmat4M(1, 1, 3, [.299, .587, .114], ROI, mask, 1);
                    [~, dat, ~, ~, ~, R{a,2}, ~, ~, ~, R{a,6}, R{a,7}, ~] = ...
                        sfrmat5M(1, 1, ROI, mask,5);

                    if ~isempty(dat)
                        if ~isnan(dat(:,end))
                            % Radial Distance: 
                            % [xleft, yTop, width, hight, Radial Distance]
    %                         BB = r{a,4};
                          
                            R{a,3}= RadialDist(r{a,4}, [size(i,2), size(i,1)]);
        
                            % ESF FWHM
                            R{a,4}=LSFfwhm(R{a,2});

%                             % Contrast
%                             R{a,5} = r{a,5};

                            if raw==0
                                mq=interp1(dat(:,1), dat(:,end), uq, 'pchip');
                            else
                                % Half the Frquency - pixels 2x apart
                                mq=interp1(dat(:,1)/2, dat(:,end), uq, 'pchip');
                            end
                            % Fitting 4th Order Polynomial
                            p = polyfit(uq,mq,4);
                            f = polyval(p,uq);
                            %Fitting Error
                            FE=abs(mq-f);
                            BadMTF=find(FE>0.1); % Fitting Thresh
                            emp2=isempty(BadMTF);
                            if emp2==1 || size(BadMTF,1)<=5 % 5=10%
                                Dat = [uq,mq]; 
                                R{a,1}=Dat;
                            else
                                R{a,1}=[];
                                R{a,2}=[];
                                R{a,3}=[];
                                R{a,4}=[];
                                R{a,5}=[];
                                R{a,6}=[];
                                R{a,7}=[];
                                R{a,8}=[];
    
                            end 
                        else 
                                R{a,1}=[];
                                R{a,2}=[];
                                R{a,3}=[];
                                R{a,4}=[];
                                R{a,5}=[];
                                R{a,6}=[];
                                R{a,7}=[];
                                R{a,8}=[];
                        end  
                    else
                        R{a,1}=[];
                        R{a,2}=[];
                        R{a,3}=[];
                        R{a,4}=[];
                        R{a,5}=[];
                        R{a,6}=[];
                        R{a,7}=[];
                        R{a,8}=[];
                    end
                end  
            end  
            switch O
                case 1 
                    NSSFRh = R;
                case 2
                    NSSFRv = R;
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%% ----- STORE DATA ----- %%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        Text = num2str(imfiles(A).name);
        MTF_Results{B,1} = Text;
        MTF_Results{B,2} = Im_lin;
        MTF_Results{B,3} = ImageInfo;
        MTF_Results{B,4} = NSSFRh;
        MTF_Results{B,5} = NSSFRv;
    end
    % Save mtf_Results cell array as .mat file in 'Results' Folder
    AA=num2str(A+ContinueA);
    filename = [resultdir '/Image-' AA '.mat'];
    parsave(filename, MTF_Results)

    time(A, 1)=toc;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% List of image names
ex = exist('namesIndex', 'var');
if ex==0
    namesIndex=zeros(imnumber, 1);
    namesIndex=num2cell(namesIndex);
end
for A=(1+ContinueA):(imnumber+ContinueA)
    AA=num2str(A);
    namesIndex{A,1}=['Image-' AA '.mat'];
end 
save([resultdir '/ImageNamesIndex.mat'], 'namesIndex', 'time')
disp('Compleated NS-SFR Extraction');
% clearvars
beep