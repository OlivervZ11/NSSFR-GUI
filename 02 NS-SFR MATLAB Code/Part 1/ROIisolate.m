function ROIs=ROIisolate(Im, BW, Ma, SE, Con, MaxEdge)
% ROIisolate is frunction to isolate the crops ising the edge isolation
% binary image. Binary Boxing is used to for this purpopse. 
% 
% Input:
%       Im:     Orginal (liner) image
%       BW:     The Edge Location Binary image output (Horizontal OR
%               Vertical only)
%       M:      Edge Mask
%       SE:    step edges - +1 = positive SE, -1 = negative SE
%       Con:    (Optional) Contrast limits [min, max], values should be 
%               between 0-1 (defult is [0, 1]).
%   MaxEdge:    (Optional) Maximum edge height - defult 128 pixels
% Outout:
%       ROIs:   A cell array containing isolated ROIs.
%
% Copyright (c) 2022 O. van Zwanenberg
% UNIVERSITY OF WESTMINSTER 
%              - COMPUTATIONAL VISION AND IMAGING TECHNOLOGY RESEARCH GROUP

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% ----- Devide large ROIs ----- %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('MaxEdge','var')
    MaxEdge = 128;
end

% Boundry Box (ROIs):
bb = regionprops(BW,'BoundingBox');
BB=struct2cell(bb);

BB=BB';

K2=0;
for k = 1 : length(bb)
    thisBB=bb(k).BoundingBox;
    if thisBB(1,4)>MaxEdge
        K1=k+K2;
        % Extract the BoundingBox coordinats
        n=thisBB(1,1);
        m=thisBB(1,2);
%             x=thisBB(1,3);
        y=thisBB(1,4);
        
        % Edge of Intrest:
        roi = imcrop(BW, thisBB);
% % % %             mask= imcrop(M{1,B},  thisBB);
        % Remove any 'non-main' edges in ROI
        lr = bwlabel(roi);
        rois=regionprops(roi,'PixelIdxList');
% % % %             lm = bwlabel(roi);
% % % %             masks=regionprops(mask,'PixelIdxList');
        % Sort based on number of pixels in each connected component:
        d = cellfun('length',{rois(:).PixelIdxList}); 
        [~,order] = sort(d,'descend');
        % Select only longest edge
        canny = ismember(lr,order(1,1)); 
        yend=find(canny(end,:)==1, 1);
        if isempty(yend)
            y=y-1;
        end
        a=y/MaxEdge;
        Floor=floor(a);
        roinum=ceil(a);
        Overlap=a-Floor;
        Overlap=1-Overlap;
        Overlap=Overlap/Floor;
        Overlap=Overlap*MaxEdge;
        
        NewBBs=zeros(roinum,1);
        NewBBs=num2cell(NewBBs);
        
        for ROInum=1:roinum
            if ROInum==1
                M=0;
                N1=find(canny(1,:)==1);
                emp=isempty(N1);
                if emp==1
                    nn=2;
                    while emp==1
                        N1=find(canny(nn,:)==1);
                        nn=nn+1;
                        emp=isempty(N1);
                    end
                end
                N1=N1(1,1);
                N=N1;
                Y=MaxEdge;
                X1=find(canny(round(Y),:)==1);
                X=abs(N-X1(1,1));
                
                if N1>X1(1,1)
                    newbb=[(N+n)-X, M+m, X, Y];
                else
                    newbb=[N+n, M+m, X, Y];
                end
                
                NewBBs{1,1}=newbb; 
            else
                M=(M+Y)-Overlap;
                N1=find(canny(round(M),:)==1);
                if emp==1
                    nn=2;
                    while emp==1
                        N1=find(canny(nn,:)==1);
                        nn=nn+1;
                        emp=isempty(N1);
                    end
                end
                N=N1(1,1);
                Y=MaxEdge;
                X1=find(canny(round(M+Y),:)==1);
                X=abs(N-X1(1,1));
                
                if N>X1(1,1)
                    newbb=[(N+n)-X, M+m, X, Y];
                else
                    newbb=[N+n, M+m, X, Y];
                end
                NewBBs{ROInum,1}=newbb; 
            end   
        end
        
        % Incert the new ROI into BB cell array, removing the large ROI
        % Incert using insertrows (Copyright (c) 2016, Jos van der Geest)
        BB = insertrows(BB, NewBBs, K1);
        BB{K1} = [];
        K2=K2+length(NewBBs);            
    end
end
%--------------------------------------------------------------------------
% Create Cell array for ROI data
l =length(BB);
ROIs=cell(l,5);
K = 0;
for k = 1 : length(BB)
    thisBB=BB{k};
    emp=isempty(thisBB);
    if emp~=1
        % Add 10 pixels either side of the edge to all
        thisBB(1,1)=thisBB(1,1)-10;
        thisBB(1,3)=thisBB(1,3)+20;

        % Minimum hight = 20 pixels
        if thisBB(1,4)<20 
            continue 
        end
        %------------------------------------------------------------------
        % Edge Location Binary Image ROI:
        canny=imcrop(BW, thisBB);
        % Remove empty rows 
        yend=find(canny(end,:)==1, 1);
        if isempty(yend)
            thisBB(1,4)=thisBB(1,4)-1;
            canny=imcrop(BW, thisBB);
        end
        % Remove any 'non-main' edges in ROI 
        l = bwlabel(canny);
        rois=regionprops(canny,'PixelIdxList');

        % Sort based on number of pixels in each connected component:
        d = cellfun('length',{rois(:).PixelIdxList}); %max number of diaginal pixels in each region
        [~,order] = sort(d,'descend');

        % Select only longest edge
        cannyS = ismember(l,order(1,1)); 
% % % %         % If there are multible same length and the first one is on the edge
% % % %         % Move to next longest edge
% % % %         [~,x]=find(cannyS==1);
% % % %         X=find(x~=1, 1);
% % % %         emp=isempty(X);
% % % %         if emp==1
% % % %             if size(order,2)>=2
% % % %                 cannyS = ismember(l,order(1,2)); 
% % % %             else
% % % %                 continue
% % % %             end
% % % %         end

% % % %         % Ensure there is a edge througout the ROI
% % % %         S = sum(cannyS,2);
% % % %         if S(1,:)==0 % Remove first
% % % %             thisBB(1,2)=thisBB(1,2)+1;
% % % %             thisBB(1,4)=thisBB(1,4)-1;
% % % %             % Resize CannyS
% % % %             cannyS = cannyS(2:size(cannyS,1),1:size(cannyS,2));
% % % %         elseif S(end,:)==0 % Remove last
% % % %             thisBB(1,4)=thisBB(1,4)-1;
% % % %             % Resize CannyS
% % % %             cannyS = cannyS(1:size(cannyS,1)-1,1:size(cannyS,2));
% % % %         elseif sum(S)~=size(cannyS,1)
% % % %              continue
% % % %         end
        %------------------------------------------------------------------
        % is Step-Edge? - 0 no SE, 1 pos SE, -1 neg SE
        % pos/neg cancel out - need more than 50% of ROI to have dominat SE
        ed=imcrop(SE, thisBB);
        ed = ed.*cannyS;
        stepedgetotal = sum(sum(ed));
        % Detemine the minimum number of confermed Step-Edges within the ROI
        SEnum=0.5*thisBB(1,4);
        %------------------------------------------------------------------
        if abs(stepedgetotal)>=SEnum % yes = step-edge
            % Image ROI:
            roi=imcrop(Im, thisBB);
            %--------------------------------------------------------------
            % Edge Mask
            mask=imcrop(Ma, thisBB);
            % Remove any 'non-main' edges in ROI 
            l = bwlabel(mask);
            masks=regionprops(logical(mask),'area');
    
            % Sort based on number of pixels in each connected component:
            d = [masks.Area]; %max number of diaginal pixels in each region
            [~,order] = sort(d,'descend');
    
            % Select only longest edge
            maskS = ismember(l,order(1,1)); 
            %--------------------------------------------------------------
            % Measure and limit contrast
            roiArea= im2double(im2gray(roi)).*double(mask);
            p95 = prctile(nonzeros(roiArea),95); 
            p05 = prctile(nonzeros(roiArea), 5);
            con = (p95-p05)/(p95+p05);
            if Con(1,1)~=0 || Con(1,2)~=1
                if con>=Con(1,2) || con<=Con(1,1) 
                    continue
                end
            end
            %--------------------------------------------------------------
            % Save values
            K=K+1;
            ROIs{K,1}=roi; % ROI
            ROIs{K,2}=cannyS; % Edge location BW
            ROIs{K,3}=maskS; % Edge Mask
            ROIs{K,4}=thisBB; % ROI Location
            ROIs{K,5}=con; % ROI Location
        else
            continue
        end
    end
end