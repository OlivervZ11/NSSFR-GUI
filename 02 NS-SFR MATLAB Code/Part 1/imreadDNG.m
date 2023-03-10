function [Ir, Ig1, Ig2, Ib, RGB]=imreadDNG(dng, display)  
    % Reads DNG file converted for .NEF file (Nikon D800 RAW file), other
    % Camera manufacturers and other Nikon models may need this code to be
    % adapted as EXIF data will be different. The mosaic sensor alignment  
    % is set to be 'rggb'.
    %
    % Code based on Steve Eddins Mathworks post, March 8, 2011.
    % Available at: https://blogs.mathworks.com/steve/2011/03/08/tips-for-reading-a-camera-raw-file-into-matlab/
    %
    % Input:
    %   dng              =  The .DNG file converted form a Nikon D800 .NEF 
    %                       file
    %   display          =  Disply [Ir, Ig1, Ig2, Ib, RGB]? 1 = yes, 0 = no
    %                       (Optional), defult = 0
    %
    % Output:
    %   Ir               =  The Red Chanel Colour Filter Array
    %   Ig1              =  The Green Chanel Colour Filter Array in row
    %                       with the Red CFA
    %   Ig2              =  The Green Chanel Colour Filter Array in row
    %                       with the Blue CFA
    %   Ib               =  The Blue Chanel Colour Filter Array
    %   RGB              =  The demosaiced image and with gamma applied
    
    ex = exist('display', 'var');
    if ex == 0 
        display = 0;
    end

    SensorAlignment = 'rggb';

    T=Tiff(dng,'r'); 
    info = imfinfo(dng);
    offsets = getTag(T,'SubIFD');
    setSubDirectory(T,offsets(1));
    % This is the Bayer pattern picture (mosaic'd)
    cfa14 = read(T); 
    cfa14 = double (cfa14);
    % Convert from 14-bit values to 16-bit values
    %     cfa=cfa14./(2^14);
        % Using taged satuation level
        WL=info.SubIFDs{1, 1}.WhiteLevel;
        cfa=cfa14./WL;
    cfa16=cfa.*(2^16);
    cfa16 = uint16(cfa16);

%     figure, imshow(RGB)
    % Seperate the RGB Channels    
[m,n] = size(cfa14);

            A1=1:2:m; % odds
            A2=2:2:m; % evens
            B1=1:2:n; % odds
            B2=2:2:n; % evens
    
    Ir  = cfa16(A1,B1);
    Ig1 = cfa16(A2,B1);
    Ig2 = cfa16(A1,B2);
    Ib  = cfa16(A2,B2);

    % White balence RGB output - each channel the same mean value
    RGB = demosaic(cfa16,SensorAlignment);
    grayImage    = double(rgb2gray(RGB));
    redChannel   = double(RGB(:, :, 1));
    greenChannel = double(RGB(:, :, 2));
    blueChannel  = double(RGB(:, :, 3));
    
    meanR = mean2(redChannel);
    meanG = mean2(greenChannel);
    meanB = mean2(blueChannel);
    meanGray = mean2(grayImage);

    redChannel = uint16(redChannel * meanGray / meanR);
    greenChannel = uint16(greenChannel * meanGray / meanG);
    blueChannel = uint16(blueChannel * meanGray / meanB);
    
    RGB = cat(3, redChannel, greenChannel, blueChannel);

    RGB=lin2rgb(RGB);
    % Add contrast 
% % % % %     RGB = histeq(RGB); 
        
    if display ==1
        figure, imshow(RGB)
        figure, imshow(Ir,[])
        figure, imshow(Ig1,[])
        figure, imshow(Ig2,[])
        figure, imshow(Ib,[])
    end
end