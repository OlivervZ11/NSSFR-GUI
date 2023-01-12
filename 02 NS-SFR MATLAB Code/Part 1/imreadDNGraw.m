function [Ir, Ig1, Ig2, Ib]=imreadDNGraw(SensorAlignment, dng)  
    % Reads DNG file converted for .NEF file (Nikon D800 RAW file),
    % other Camera manufacturers and other Nikon models may need 
    % this code to be adapted as EXIF data will be different.
    %
    % Code based on Steve Eddins Mathworks post, March 8, 2011.
    % Available at: https://blogs.mathworks.com/steve/2011/03/08/
    % tips-for-reading-a-camera-raw-file-into-matlab/
    %
    % Input:
    %   dng              =  The .DNG file (Tested via converted   
    %                       Nikon D800 .NEF files)
    %
    % Output:
    %   Ir               =  The Red Chanel Colour Filter Array
    %   Ig1              =  The Green Chanel Colour Filter Array 
    %                       in row with the Red CFA
    %   Ig2              =  The Green Chanel Colour Filter Array 
    %                       in row with the Blue CFA
    %   Ib               =  The Blue Chanel Colour Filter Array
    
%     SensorAlignment = app.BayerPatternDropDown.Value;

    T=Tiff(dng,'r'); 
    info = imfinfo(dng);
    offsets = getTag(T,'SubIFD');
    setSubDirectory(T,offsets(1));
    % Convert from 14-bit values to 16-bit values
    % This is the Bayer pattern picture (mosaiced)
    cfa14 = read(T); 
    cfa14 = double (cfa14);

    % Using taged satuation level
    WL=info.SubIFDs{1, 1}.WhiteLevel;
    cfa=cfa14./WL;
    cfa16=cfa.*(2^16);
    cfa16 = uint16(cfa16);

    % Seperate the RGB Channels    
    [m,n] = size(cfa14);

    A1=1:2:m; % odds
    A2=2:2:m; % evens
    B1=1:2:n; % odds
    B2=2:2:n; % evens
    
    switch SensorAlignment
        case 'rggb'
            Ir  = cfa16(A1,B1);
            Ig1 = cfa16(A2,B1);
            Ig2 = cfa16(A1,B2);
            Ib  = cfa16(A2,B2);
        case 'gbrg'
            Ig1  = cfa16(A1,B1);
            Ib = cfa16(A2,B1);
            Ir = cfa16(A1,B2);
            Ig2  = cfa16(A2,B2);
        case 'grbg'
            Ig1  = cfa16(A1,B1);
            Ir = cfa16(A2,B1);
            Ib = cfa16(A1,B2);
            Ig2  = cfa16(A2,B2);
        case 'bggr'
            Ib  = cfa16(A1,B1);
            Ig1 = cfa16(A2,B1);
            Ig2 = cfa16(A1,B2);
            Ir  = cfa16(A2,B2);
    end
end 