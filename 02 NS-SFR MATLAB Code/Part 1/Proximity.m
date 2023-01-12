function [Proximity]=Proximity(BW, Coordinates, lowerLimit, upperLimit)

% PROXIMITY A proximity filter that takes either the coordinates from the
% vertical or horizontal edges, then measures and limits the proximity of 
% the detected edges.
% Any edges between the proximity of lowerLimit and upperLimit are removed 
% from the Vertical/Horizontal Coordinate list. 
%
% Input: 
%    BW:            Initial Binary Image
%    Coordinates:   List of Vertical/Horizontal Coordinates
%    lowerLimit:    The lower threshold, any edges above this threshold are
%                   eliminated (as long as the proximity does not exceed 
%                   the upperLimit).
%    upperLimit:    The upper threshold, any edges bellow this threshold 
%                   are eliminated (as long as as the proximity is above  
%                   that of the lowerLimit).   
%
% Output:
%    Proximity:                 A new binary image, where the edges that 
%                               are within a set proximity are removed.
% 
% Copyright (c) 2022 O. van Zwanenberg, University of Westminster


%--------------------------------------------------------------------------
Co=zeros(size(Coordinates,1),3);
Co(:,1:2)=sortrows(Coordinates,2); 
% allIdx = (1:(size(Coordinates,1)))';
Max=max(Coordinates(:,2));
Proximity=zeros(size(BW(:,:)));
for C=1:Max
     [R, ~] = find(Co(:,2) == C);
%     R=(Co(:,2) == C);
%     R = allIdx(R ~= 0) ;
%     if isempty (R)
%         R=[];
%     end

    % If one coordinate in row
    if size(R,1)==1
        Co(R,3) = 1;
        continue
    end
    CUT = zeros(size(R,1), 3);
    CUT(:,1:2)=Co(R, 1:2);
 
    %Determine if HCECxcut is empty - 1=true, 0=false
    if isempty(CUT)==0
        Co2=CUT(2:end, 1);
        proximity = Co2-Co(R(1:end-1),1);
        P = find(proximity>=upperLimit | proximity<=lowerLimit);
        if isempty(P)==1
            continue
        end
        CUT(P, 3)=1;
        if P(end)==size(P,1)
            CUT(end, 3)=1;
        end
        Z = find(CUT(:, 3)==0);
        if isempty(Z)==0
            CUT(Z+1, 3)=0;
        end
        R=R(CUT(:, 3)==1);
        Co(R,3) = 1;

%         for E=1:size(CUT,1)-1
%             x1=CUT(E,1);
%             G=E+1;
%             x2=CUT(G,1);      
%             proximity=(x2-x1);
%             
%             %create a new Horizoantal Canny Edge Detection 
%             %- with proximity reduction
%             if proximity>=upperLimit || proximity<=lowerLimit
%                 Proximity(C, x1)=1;
%                 Proximity(C, x2)=1;
%             else
%                 Proximity(C, x1)=0;
%                 Proximity(C, x2)=0;
%             end
%         end
    end
end
Co = Co(any(Co(:,3),2),:);
Proximity(sub2ind(size(Proximity),Co(:,2),Co(:,1))) = 1;
Proximity=logical(Proximity);