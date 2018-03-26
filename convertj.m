function [n, m, sc] = convertj(j, ztype)
%CONVERTJ Convert ordering parameter j to n, m, sc parameters.
%
% function [n, m, sc] = convertj(j, ztype)
%
%   INPUT:
%     j     = Zernike polynomial order number.  j >= 0 and integer.
%             For 'FRINGE' 1 <= j <= 37.
%             For 'ISO' 0 <= j <= 35.
%     ztype = type of odering: 'FRINGE', 'ISO', 'WYANT', 'B&W',
%                              'STANDARD', 'MAHAJAN', 'NOLL',
%                              'HEXAGON', 'HEXAGONR30', 'ELLIPSE',
%                              'RECTANGLE', 'SQUARE', 'ANNULUS'
%
%   OUTPUT:
%     n  = radial order.
%     m  = azimuthal order
%     sc = indicate whether sine or cosine or neither factor
%          's' means sin() factor
%          'c' means cos() factor
%          ' ' means m = 0 so has no sin() or cos() factor
%

%
% DANGER: Validation of input is NOT performed.
%         This function should only be called from ZernikeCalc
%         which does the input validation.
%


switch ztype
    case 'FRINGE'
        [n,m,sc] = fringe(j);
        
    case {'ISO', 'WYANT'}
        [n,m,sc] = fringe(j+1);
        
    case {'B&W', 'STANDARD'}
        [n,m,sc] = bw(j);
        
    case {'MAHAJAN','NOLL', ...
            'HEXAGON', 'HEXAGONR30', 'ELLIPSE', ...
            'RECTANGLE', 'SQUARE', 'ANNULUS'}
        [n,m,sc] = mahajan(j);
        
end % switch ztype

end % function convertj