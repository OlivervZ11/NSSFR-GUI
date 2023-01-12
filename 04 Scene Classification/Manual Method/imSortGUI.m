% GUI to go through folder of images and lable their content manually
% Copyright (c) 2023 O. van Zwanenberg
% UNIVERSITY OF W1ESTMINSTER 
%              - COMPUTATIONAL VISION AND IMAGING TECHNOLOGY RESEARCH GROUP

clc; close all; clear all;

% USer input: Continue?
answer = questdlg('Would you like to build on a previous imSort?', ...
	'Continue?', ...
	'Start a New Dataset','Continue and load previous imSort','Start a New Dataset');
% Handle response
switch answer
    case 'Start a New Dataset'
        SceneClass=[];
        answer2 = questdlg('Would you like to upload Dropdown options?', ...
        'Upload Dropdown?', ...
	    'No','Yes','No');
        % Handle response
        switch answer2
             case 'No'
                 Dropdowns=[];
            case 'Yes'
               [file,path] = uigetfile('*.mat');
               load([path file]);
        end 
    case 'Continue and load previous imSort'
        %User select folder of images
        ConPath = uigetdir;
        load([ConPath '/SceneClassData.mat']);
        load([ConPath '/Dropdowns.mat']);
end

% RAW or TIFF
answer = questdlg('Are te images TIFF or DNG (RAW) files?', ...
	'File Type', ...
	'TIFF','DNG','TIFF');
% Handle response
switch answer
    case 'TIFF'
        RAW=0;
    case 'DNG'
        RAW=1;

end

%User select folder of images
path = uigetdir;
if RAW == 0
    filePattern = fullfile(path, '*.tif');
else
    filePattern = fullfile(path, '*.dng');
end
ImFiles = dir(filePattern);

names={ImFiles.name};
imSort(path, names, SceneClass, Dropdowns, RAW)


% Create UIFigure and components
function imSort(path, names, SceneClass, Dropdowns, RAW)
    % Obtain image data
    ImNum=size(names,2);
    k=1;
    i = fullfile(path,names{k});
    if RAW==0
        im=imread(i);
    else
        [~, ~, ~, ~, im]=imreadDNG(i, 0);
        
    end
    
    MD=imfinfo(i);
    iso=MD.DigitalCamera.ISOSpeedRatings;

    % Create UIFigure and hide until all components are created
    UIFigure = uifigure('Visible', 'off');
    UIFigure.Position = [100 100 644 423];
    UIFigure.Name = 'imSort';
    
    % Create ImageAxes
    ImageAxes = uiaxes(UIFigure);
    ImageAxes.XTick = [];
    ImageAxes.XTickLabel = {'[ ]'};
    ImageAxes.YTick = [];
    ImageAxes.Position = [15 85 431 305];
    Im=image(im,'Parent', ImageAxes);
    set(ImageAxes,'visible','off');

% % %     % Create PREVIOUSButton
% % %     PREVIOUSButton = uibutton(UIFigure, 'push');
% % %     PREVIOUSButton.Position = [95 25 115 37];
% % %     PREVIOUSButton.Text = 'PREVIOUS';
% % %     PREVIOUSButton.ButtonPushedFcn = @changeImNeg;

    % Create NEXTButton
    NEXTButton = uibutton(UIFigure, 'push');
    NEXTButton.Position = [236 25 115 37];
    NEXTButton.Text = 'NEXT';
    NEXTButton.ButtonPushedFcn = @changeImPos;
    
    % ---------------------------------------------------------------------
    % Lighting:
    
    % Create DropDownLabel
    DropDownLabel = uilabel(UIFigure);
    DropDownLabel.HorizontalAlignment = 'right';
    DropDownLabel.Position = [449 354 66 22];
    DropDownLabel.Text = 'Drop Down';

    % Create DropDown
    DropDown = uidropdown(UIFigure);
    DropDown.Position = [530 354 100 22];
    DropDown.Items = {'Other'};
    if isempty(Dropdowns)
         DropDown.Items = {'Other'};
    else
         DropDown.Items = (['Other'; Dropdowns{1,2}(:,1)]);
    end


    % Create OtherEditFieldLabel
    OtherEditFieldLabel = uilabel(UIFigure);
    OtherEditFieldLabel.HorizontalAlignment = 'right';
    OtherEditFieldLabel.Position = [476 321 39 22];
    OtherEditFieldLabel.Text = 'Other:';

    % Create OtherEditField
    OtherEditField = uieditfield(UIFigure, 'text');
    OtherEditField.Position = [530 321 100 22];
    
    % ---------------------------------------------------------------------
    % Location:
    % Create DropDown_2Label
    DropDown_2Label = uilabel(UIFigure);
    DropDown_2Label.HorizontalAlignment = 'right';
    DropDown_2Label.Position = [449 260 66 22];
    DropDown_2Label.Text = 'Drop Down';

    % Create DropDown_2
    DropDown_2 = uidropdown(UIFigure);
    DropDown_2.Position = [530 260 100 22];
   
    if isempty(Dropdowns)
         DropDown_2.Items = {'Other'};
    else
         DropDown_2.Items = (['Other'; Dropdowns{1,1}(:,1)]);
    end
    
    % Create OtherEditField_2Label
    OtherEditField_2Label = uilabel(UIFigure);
    OtherEditField_2Label.HorizontalAlignment = 'right';
    OtherEditField_2Label.Position = [476 227 39 22];
    OtherEditField_2Label.Text = 'Other:';

    % Create OtherEditField_2
    OtherEditField_2 = uieditfield(UIFigure, 'text');
    OtherEditField_2.Position = [530 227 100 22];

    % ---------------------------------------------------------------------
    % Create MainLightingConditionLabel
    MainLightingConditionLabel = uilabel(UIFigure);
    MainLightingConditionLabel.FontSize = 14;
    MainLightingConditionLabel.FontWeight = 'bold';
    MainLightingConditionLabel.Position = [449 385 172 22];
    MainLightingConditionLabel.Text = 'Main Light Source';

    % Create ImageLocationLabel
    ImageLocationLabel = uilabel(UIFigure);
    ImageLocationLabel.HorizontalAlignment = 'center';
    ImageLocationLabel.FontSize = 14;
    ImageLocationLabel.FontWeight = 'bold';
    ImageLocationLabel.Position = [449 296 172 22];
    ImageLocationLabel.Text = 'Primary Image Subject';
    % ---------------------------------------------------------------------
    % Create MainLocationButtonGroup
    MainLocationButtonGroup = uibuttongroup(UIFigure);
    MainLocationButtonGroup.Title = 'General Environment';
    MainLocationButtonGroup.Position = [500 120 123 98];

    % Create ManMadeButton
    ManMadeButton = uitogglebutton(MainLocationButtonGroup);
    ManMadeButton.Text = 'Man Made';
    ManMadeButton.Position = [11 44 100 23];
    ManMadeButton.Value = true;

    % Create IndoorButton
    IndoorButton = uitogglebutton(MainLocationButtonGroup);
    IndoorButton.Text = 'Indoor';
    IndoorButton.Position = [11 23 100 23];

    % Create NatureButton
    NatureButton = uitogglebutton(MainLocationButtonGroup);
    NatureButton.Text = 'Nature';
    NatureButton.Position = [11 2 100 23];
    
    % ---------------------------------------------------------------------
    % Create Gauge
    Gauge = uigauge(UIFigure, 'ninetydegree');
    Gauge.Position = [507 25 92 92];
    Gauge.Limits = [1 ImNum];
    Gauge.Value = 1;
    Gauge.MajorTicks = [1 ImNum];
    
    % Show the figure after all components are created
    UIFigure.Visible = 'on';

    if isempty(SceneClass)
        % Ims
        ImList = cell(size(names,2), 5);
        K3 = 0;
        % Lighting 
        K1=0;
        Li=num2cell([0,0,0,0]);
        % Location
        K2=0;
        Loc=num2cell([0,0,0,0]);
        if ~isempty(Dropdowns)
            % Location
            K2 = size(Dropdowns{1, 1},1);
            for d = 1:K2
                Loc{d,1}=Dropdowns{1, 1}{d,1};
                Loc{d,2}=0;
                Loc{d,3}=0;
                Loc{d,4}=0;
            end
            % Lighting 
            K1 = size(Dropdowns{1, 2},1);
            for d = 1:K1
                Li{d,1}=Dropdowns{1, 2}{d,1};
                Li{d,2}=0;
                Li{d,3}=0;
                Li{d,4}=0;
            end
        end
        
    else
        % Lighting 
        Li = SceneClass{1,3};
        K1 = size(SceneClass{1,3},1);
        % Location
        Loc= SceneClass{1,2};
        K2 = size(SceneClass{1,2},1);
        % Ims
        ImList = SceneClass{1,1};
        K3 = size(SceneClass{1,1},1);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%% ----- FUNCTIONS ----- %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function changeImPos(~,~,~)

        % Save Lables:

        %Manmade, Indoor or Nature
        if ManMadeButton.Value == true
             MNI = 2;
        elseif IndoorButton.Value == true
             MNI = 3;
        elseif NatureButton.Value == true
             MNI = 4;
        end

        % Lighting Array
        % New Lable
        DropLable=get(DropDown, 'Value');
        if strcmp(DropLable,'Other')
            NewLable=get(OtherEditField, 'Value');
             set(OtherEditField, 'Value', '');
            DropLable = NewLable;
            Index = find(strcmp(Li, NewLable));
            emp=isempty(Index);
            % Check to see if submission is New
            if emp==1
                K1=K1+1;
                Li{K1,1}=NewLable;
                Li{K1,2}=0;
                Li{K1,3}=0;
                Li{K1,4}=0;
                Li{K1,MNI}=1;
                Droplist=get(DropDown, 'Items');
                numList=size(Droplist,2);
                Droplist{1,numList+1}=NewLable;
                % sort alphabetically
                [~, ix] = sort(Droplist(1, 2:end));
                Droplist(:, 2:end) = Droplist(:, ix+1);
                set(DropDown, 'Items', Droplist);
            else
                % Add to a Lable
                Li{Index,MNI}=Li{Index,MNI}+1;
            end
        else
            % Add to a Lable
            Index = find(strcmp(Li, DropLable));
            Li{Index,MNI}=Li{Index,MNI}+1;

        end
        li = DropLable;

        % Location 
        % New Lable
        DropLable=get(DropDown_2, 'Value');
        if strcmp(DropLable,'Other')
            NewLable=get(OtherEditField_2, 'Value');
            set(OtherEditField_2, 'Value', '');
            DropLable = NewLable;
            Index = find(strcmp(Loc, NewLable));
            emp=isempty(Index);
            % Check to see if submission is New
            if emp==1
                K2=K2+1;
                Loc{K2,1}=NewLable;
                Loc{K2,2}=0;
                Loc{K2,3}=0;
                Loc{K2,4}=0;
                Loc{K2,MNI}=1;
                Droplist=get(DropDown_2, 'Items');
                numList=size(Droplist,2);
                Droplist{1,numList+1}=NewLable;
                % sort alphabetically
                [~, ix] = sort(Droplist(1, 2:end));
                Droplist(:, 2:end) = Droplist(:, ix+1);
                set(DropDown_2, 'Items', Droplist);
            else
                % Add to a Old Lable
                Loc{Index,MNI}=Loc{Index,MNI}+1;
            end
        else
            % Add to a Old Lable
            Index = find(strcmp(Loc, DropLable));
            Loc{Index,MNI}=Loc{Index,MNI}+1;
        end
        Lo = DropLable;

        % Image list
        ImList{k+K3,1} = names{k};
        ImList{k+K3,2} = Lo;
        ImList{k+K3,3} = li;
        switch MNI
            case 2
                ImList{k+K3,4} = 'ManMade';
            case 3
                ImList{k+K3,4} = 'Indoor';
            case 4
                ImList{k+K3,4} = 'Nature';
        end
        ImList{k+K3,5} = iso;
            
        % Display next image
        k=k+1;
        if k<=ImNum
            i = fullfile(path,names{k});
            if RAW==0
                im=imread(i);
            else
                [~, ~, ~, ~, im]=imreadDNG(i, 0);
                MD=imfinfo(i);
                iso=MD.DigitalCamera.ISOSpeedRatings;
            end
            set(Im,'CData',im);
            Gauge.Value = k;
        else
%             k=k-1;
            % End of all images - Save results
            SceneClass=cell(1,3);
            SceneClass{1,1}=ImList;
            SceneClass{1,2}=Loc;
            SceneClass{1,3}=Li;
            % Delete Fig ???
%             close (gcf)
            % Save
            savepath = uigetdir;
            filename = [savepath '/SceneClassData.mat'];
            save(filename, 'SceneClass')
%             uisave({'SceneClass'},'SceneClassData')

            % ImList Table
            filename = [savepath '/SceneClass.xlsx'];
            
            Filenames = ImList(:,1);
            Location = ImList(:,2);
            Lighting = ImList(:,3);
            Class =  ImList(:,4);
            ReportedISO =  ImList(:,5);
            ImListTable = table(Filenames,Location,Lighting,Class, ReportedISO);
            writetable(ImListTable,filename,'Sheet',1)
           
            Subject = Loc(:,1);
            ManMade = Loc(:,2);
            Indoor =  Loc(:,3);
            Nature =  Loc(:,4);
            LocListTable = table(Subject,ManMade,Indoor,Nature);
            writetable(LocListTable,filename,'Sheet',2)

            Lighting = Li(:,1);
            ManMade = Li(:,2);
            Indoor =  Li(:,3);
            Nature =  Li(:,4);
            LiListTable = table(Lighting,ManMade,Indoor,Nature);
            writetable(LiListTable,filename,'Sheet',3)

            % sort alphabetically
            [~, ix] = sort(Subject(1:end, 1));
            Subject(1:end,:) = Subject(ix,:);

            [~, ix] = sort(Lighting(1:end, 1));
            Lighting(1:end,:) = Lighting(ix,:);
            
            Dropdowns  =cell(1,2);
            Dropdowns{1,1} = Subject;
            Dropdowns{1,2} = Lighting;
            filename = [savepath '/Dropdowns.mat'];
            save(filename, 'Dropdowns');
        end
        if k==ImNum
            set(NEXTButton,'Text', 'SAVE');
            
        end
    end
end
% Code that executes before app deletion
function delete(UIFigure)

    % Delete UIFigure when app is deleted
    delete(UIFigure)
end