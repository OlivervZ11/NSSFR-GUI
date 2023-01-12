% Manually select ROIs from an image and measure the ISO12233:2017 e-SFR

% Copyright (c) 2023 O. van Zwanenberg
% UNIVERSITY OF WESTMINSTER 

clc; close all; clear all;

% Radial Annuli
RA = 3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Select Image
answer = questdlg('Test Chart Capture File Format', ...
	'Data', ...
	'RAW(.dng)','TIFF(.tif)', 'RAW(.dng)');
switch answer
    case 'RAW(.dng)' 
        raw = 1;
        rgb = 4; %rggb
    case 'TIFF(.tif)'
        raw = 0;
        rgb=1;
end

imgPath = uigetdir;
if raw==0
images  = dir([imgPath '/*.tif']);
else
images  = dir([imgPath '/*.dng']);
end
N = length(images);
IM{N,rgb} = [];
for idx = 1:N
    switch raw
        case 0
            IM{idx,1} = imread([imgPath '/' images(idx).name]);
        case 1
            [Ir, Ig1, Ig2, Ib, Irgb]=imreadDNG([imgPath '/' images(idx).name], 0);
            IM{idx,1} = Ir;
            IM{idx,2} = Ig1;
            IM{idx,3} = Ig2;
            IM{idx,4} = Ib;
            IM{idx,5} = imresize(Irgb,size(Ig1));
    end
end
% Im=imread([list(x).name]);

% Data Cell Arrays - format matches the Framework output
ISO12233Chart = zeros (rgb,6);
ISO12233Chart=num2cell(ISO12233Chart);

Dat = zeros (1,12);
Dat=num2cell(Dat);

for RGB = 1:rgb
ISO12233Chart{RGB,4}=Dat;
ISO12233Chart{RGB,5}=Dat;
end
% -------------------------------------------------------------------------
% Select the ROIs from Im
ISO12233Chart=SelectROIs(IM, ISO12233Chart, raw); 
% ISO12233Chart{1,1}=file;
switch raw
    case 0
        ISO12233Chart{1,2}=IM;
    case 1
        ISO12233Chart{1,2}=Ir;
        ISO12233Chart{2,2}=Ig1;
        ISO12233Chart{3,2}=Ig2;
        ISO12233Chart{4,2}=Ib;
end

% Save the ISO12233Chart.mat
% % % answer = questdlg('Would you like a save the SFR data in .xlsx?', ...
% % % 	'Save data', ...
% % % 	'Yes', 'No', 'Yes');
% % % % Handle response
% % % switch answer
% % %     case 'Yes'
% % % %         uisave('ISO12233Chart','ISO12233ChartData');
        % Transfer the ISO12233 e-SFR data in Excel Doc. 
        savepath = uigetdir(imgPath,'Save Data in .xlsx Documment');
        freq=(0:0.01:0.50)';
        
        if size(ISO12233Chart,1)==1
            raw=0;
            rgb=1;
            rggb=1;
        else
            raw=1;
            rgb=3; 
            rggb=4;
        end
        % Radial Annuli  
        im = ISO12233Chart{1, 2} ; 
        [R, ~] = RadAnnuli(im{1, 1}   , RA);
        
        eSFR = cell(rgb,RA+1);
        LSF  = cell(rgb,RA);
        
        
        
        % Sort Data acording to annuli 
        for O = 4:5 % H/V edges
            for RGGB = 1:rggb
                if RGGB<=2
                   RGB=RGGB;
                   x = zeros(1,RA+1)+1;
                else
                   RGB=RGGB-1;
                   if RGGB==4
                        x = zeros(1,RA+1)+1;
                   end
                end
                for A = 1:size(ISO12233Chart{RGGB, O},1)
                    r = ISO12233Chart{RGGB, O}{A, 10};
                    for B = 1:RA-1
                        if r<=R(1,B)
                            eSFR{RGB,B}(:,x(B)) = ISO12233Chart{RGGB, O}{A, 11}{1, 1}(:,2);
                            LSF{RGB,B}{x(B),2}  = ISO12233Chart{RGGB, O}{A, 11}{1, 2};
                            % Add all eSFRs to the entire frame cell
                            eSFR{RGB,RA+1}(:,x(RA+1)) = ISO12233Chart{RGGB, O}{A, 11}{1, 1}(:,2);
        
                            x(B)=x(B)+1;
                            x(RA+1)=x(RA+1)+1;
                            continue
                        end
                    end
                    if r>R(1,RA-1)
                        eSFR{RGB,RA}(:,x(RA)) = ISO12233Chart{RGGB, O}{A, 11}{1, 1}(:,2);
                        LSF{RGB,RA}{x(RA),2}  = ISO12233Chart{RGGB, O}{A, 11}{1, 2};
                        % Add all eSFRs to the entire frame cell
                        eSFR{RGB,RA+1}(:,x(RA+1)) = ISO12233Chart{RGGB, O}{A, 11}{1, 1}(:,2);
                        x(RA)=x(RA)+1;
                        x(RA+1)=x(RA+1)+1;
                    end
                end
            end
        [aveSFR, ~]=medianSFR(LSF, [1, 0.75, 0.5], raw);
            for RGB = 1:rgb
                switch O
                    case 4
                        if raw ==0
                            filename = [savepath '/ISO12233eSFR_Vertical.xlsx'];
                        else
                            switch RGB
                                case 1
                                 filename = [savepath '/ISO12233eSFR_Vertical(R).xlsx'];
                                 case 2
                                 filename = [savepath '/ISO12233eSFR_Vertical(G).xlsx'];
                                 case 3
                                 filename = [savepath '/ISO12233eSFR_Vertical(B).xlsx'];
                            end
                        end
                    case 5
                        if raw ==0
                            filename = [savepath '/ISO12233eSFR_Horizontal.xlsx'];
                        else
                            switch RGB
                                case 1
                                 filename = [savepath '/ISO12233eSFR_Horizontal(R).xlsx'];
                                 case 2
                                 filename = [savepath '/ISO12233eSFR_Horizontal(G).xlsx'];
                                 case 3
                                 filename = [savepath '/ISO12233eSFR_Horizontal(B).xlsx'];
                            end
                        end
                end
                
                for B = 1:RA+1
                    if isempty(aveSFR{RGB, B})==1
                        continue
                    end
                    AveSFR = aveSFR{RGB, B}(:,2);
                    SFRTable = table(AveSFR, freq);
                    for A = 1:size(eSFR{1, B},2)
                        SFRTable = [ table(eSFR{1, B}(:,A), 'VariableNames', {['ROI ' num2str(A)]})  SFRTable]; 
                    end
                
                SFRTable = SFRTable(:,flip(1:end));
                writetable(SFRTable,filename,'Sheet',B)
                end
            end
        end

% %     case 'No'
% %         disp('SFR data is not saved to disk');
% % end
% -------------------------------------------------------------------------
% % % % Plot the SFR data in GUI
% % % answer = questdlg('Would you like a plot the SFR data?', ...
% % % 	'Plot data', ...
% % % 	'Yes','No');
% % % % Handle response
% % % switch answer
% % %     case 'Yes'
% % %         plotSFR(ISO12233Chart);
% % %     case 'No'
% % %         
% % % end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Selection GUI
function ISO12233Chart=SelectROIs(IM, ISO12233Chart, raw)
    % Take first image
    i=1;
    if raw == 0
        Im=IM{i,1};
        rgb=1;
    else
        Im=lin2rgb(IM{i,5}); % Irgb
        rgb=4;
    end
    imNum=size(IM,1);
    
    % ROI number indercator 
    roiNum=0;
    [y,x,~]=size(Im);
    ImS(1,1)=x;
    ImS(1,2)=y;
    
    % Orientation indercator 
    o=4;
    
%     [x,y,~] = size(Im);
    % Create UIFigure and hide until all components are created
    ROISelection = uifigure('Visible', 'off');
    ROISelection.Color = [0.651 0.651 0.651];
    ROISelection.Position = [100 100 1350 759];
    ROISelection.Name = 'ROI Selection';

    % Create Image
    Image = uiaxes(ROISelection);
    Image.Position = [59,16,862,575];
    Image.BackgroundColor = [0.651 0.651 0.651]; 
%     Image.XLim = [0 x];
%     Image.YLim = [0 y];
    set(Image,'visible','off');
    image(Im,'Parent', Image); 
    
    % Linearise the Image
    Im=rgb2lin(Im);
    
    % Create SelectButton
    SelectButton = uibutton(ROISelection, 'push');
    SelectButton.FontSize = 18;
    SelectButton.FontWeight = 'bold';
    SelectButton.Position = [1033,343,242,74];
    SelectButton.Text = 'Select';
    SelectButton.ButtonPushedFcn = @ROIselect;

    % Create NextButton
    NextButton = uibutton(ROISelection, 'push');
    NextButton.FontSize = 18;
    NextButton.FontWeight = 'bold';
    NextButton.Position = [1033 175 242 74];
    NextButton.Text = 'Next';
    NextButton.ButtonPushedFcn = @Next;
    N=1;
    
    % Create TitleLabel
    TitleLabel = uilabel(ROISelection);
    TitleLabel.HorizontalAlignment = 'center';
    TitleLabel.FontSize = 24;
    TitleLabel.Position = [59,692,862,44];
    TitleLabel.Text = 'Horizontal EDGE Selection';
        
    % Create ROILabel
    ROILabel = uilabel(ROISelection);
    ROILabel.HorizontalAlignment = 'center';
    ROILabel.FontSize = 20;
    ROILabel.Position = [59,629,862,39];
    ROILabel.Text = 'ROI Size: [100, 100, 128, 64]';

    % Create numLabel
    numLabel = uilabel(ROISelection); 
    numLabel.FontSize = 20;
    numLabel.Position = [1033,290,242,39];
    numLabel.Text = '0';
    
    % Create ROIheightSpinnerLabel
    ROIheightSpinnerLabel = uilabel(ROISelection);
    ROIheightSpinnerLabel.HorizontalAlignment = 'right';
    ROIheightSpinnerLabel.FontSize = 18;
    ROIheightSpinnerLabel.Position = [1033 591 91 23];
    ROIheightSpinnerLabel.Text = 'ROI height';

    % Create ROIheightSpinner
    ROIheightSpinner = uispinner(ROISelection);
    ROIheightSpinner.FontSize = 18;
    ROIheightSpinner.Position = [1139 590 136 24];
    ROIheightSpinner.Value = 128;
    ROIheightSpinner.ValueChangingFcn=@dispROIsize;

    % Create ROIWidthSpinnerLabel
    ROIWidthSpinnerLabel = uilabel(ROISelection);
    ROIWidthSpinnerLabel.HorizontalAlignment = 'right';
    ROIWidthSpinnerLabel.FontSize = 18;
    ROIWidthSpinnerLabel.Position = [1033 557 88 23];
    ROIWidthSpinnerLabel.Text = 'ROI Width';

    % Create ROIWidthSpinner
    ROIWidthSpinner = uispinner(ROISelection);
    ROIWidthSpinner.FontSize = 18;
    ROIWidthSpinner.Position = [1139 556 136 24];
    ROIWidthSpinner.Value = 64;
    ROIWidthSpinner.ValueChangingFcn=@dispROIsize;

    % Show the figure after all components are created
    ROISelection.Visible = 'on';
    
    rect = drawrectangle('Color',[1 0 0], 'Parent', Image, ...
        'Position', [100, 100, 128, 64]);    
    roiL = rect.Position;
    addlistener(rect,'ROIMoved',@(src, evt) roiChange(src,evt));
    
    uiwait(ROISelection);
    % ---------------------------------------------------------------------
    % GUI Functions:
    
    % ROI change via spinners
    function dispROIsize(~,~,~)
        h=get(ROIheightSpinner, 'Value');
        w=get(ROIWidthSpinner, 'Value');
        
        p=get(rect, 'Position');
        p(1,3)=h;
        p(1,4)=w;
        set(rect, 'Position', p);
        
        text=['ROI Size: [' num2str(p) ']' ];
        set (ROILabel, 'Text', text);
    end
    
    % ROI change
    function roiChange(~,~)
%         assignin('base',roi,evt.CurrentPosition);
%         roiL=evt.CurrentPosition;
        p=get(rect, 'Position');
        text=['ROI Size: [' num2str(p) ']' ];
        set (ROILabel, 'Text', text);
        
        set(ROIheightSpinner, 'Value', p(1,3));
        set(ROIWidthSpinner, 'Value', p(1,4));
    end
    
    % Select the ROI
    function ROIselect(~,~,~)
        % Mark ROI
        p=get(rect, 'Position');
        rectangle('Position', p,'Edgecolor', 'r', 'Parent', Image);
        
        roiNum=roiNum+1;
        
        % crop section
        for RGB = 1:rgb
            if raw ==0
                 ROI = imcrop(Im,p);
                 ROI = rgb2gray(ROI);
            else
                ROI = imcrop(IM{i,RGB},p);
            end
               
            % Caculate data and save
            % frequency range
            uq=(0:0.01:0.50)';

            normdist=RadialDist(p, ImS);
            [~, dat, ~, ~, lsf, ~, ~, Con, Angle, Clip] = sfrmat4(1, 1, 3, [.299, .587, .114], ROI, 1);

            if raw==0
                mq=interp1(dat(:,1), dat(:,end), uq, 'pchip');
            else
                % Half the Frquency - pixels 2x apart
                mq=interp1(dat(:,1)/2, dat(:,end), uq, 'pchip');
            end

            d(:,1)=uq;
            d(:,2)=mq;

            Data=zeros(1,2);
            Data=num2cell(Data);
            Data{1,1}=d;
            Data{1,2}=lsf;
            
            ISO12233Chart{RGB,o}{roiNum, 1}=ROI;
            ISO12233Chart{RGB,o}{roiNum, 4}=Angle;
            ISO12233Chart{RGB,o}{roiNum, 5}=Con;
            ISO12233Chart{RGB,o}{roiNum, 6}=Clip;
            ISO12233Chart{RGB,o}{roiNum, 9}=p;
            ISO12233Chart{RGB,o}{roiNum, 10}=normdist;
            ISO12233Chart{RGB,o}{roiNum, 11}=Data;
        end
        % Add to the diplayed counter
        text=num2str(roiNum);
        set (numLabel, 'Text', text);
    end

    % Move to next image or orientation and then move to finish
    function Next(~,~,~)
        i=i+1;
        if i==imNum+1
            roiNum=0;
            switch N
                case 1 % move to Vertical
                    set(TitleLabel,'Text','Vertical EDGE Selection');
                    o=5;
                    set (numLabel, 'Text', '0');
                case 2 % Finish GUI
                    uiresume (ROISelection);
                    delete(ROISelection);

            end
            if N==1
                N=N+1;
                p=get(rect, 'Position');
                i=1;
                if raw == 0
                    Im=IM{i,1};
                else
                    Im=lin2rgb(IM{i,5}); % Irgb
                end
                [y,x,~]=size(Im);
                ImS(1,1)=x;
                ImS(1,2)=y;
                image(Im,'Parent', Image); 
                % Linearise the Image
                Im=rgb2lin(Im);
                rect = drawrectangle('Color',[1 0 0], 'Parent', Image, ...
                    'Position', p);   
                addlistener(rect,'ROIMoved',@(src, evt) roiChange(src,evt));
            end
        else
            p=get(rect, 'Position');
            if raw == 0
                Im=IM{i,1};
            else
                Im=lin2rgb(IM{i,5}); % Irgb
            end
            [y,x,~]=size(Im);
            ImS(1,1)=x;
            ImS(1,2)=y;
            image(Im,'Parent', Image); 
            % Linearise the Image
            Im=rgb2lin(Im);
            rect = drawrectangle('Color',[1 0 0], 'Parent', Image, ...
                'Position', p);   
            addlistener(rect,'ROIMoved',@(src, evt) roiChange(src,evt));
            if i==imNum && o==5
               set(NextButton,'Text', 'Finish'); 
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
