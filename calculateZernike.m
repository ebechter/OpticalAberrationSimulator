function result = calculateZernike(n, m, sc, rho, theta, ZernikeDef)
%
%  Calculates the Zernike polynomial value for the given pixel location.
%
%   INPUT:
%     n           = radial order.
%     m           = azimuthal order.
%     sc          = ' ' for m = 0,  = 's' for sin() term,  = 'c' for cos() term.
%     rho         = normalized radial distance to pixel location.
%     theta       = angle from the x-axis in radians of pixel location.
%     ZernikeDef  = One of 'MAHAJAN', 'NOLL', 'FRINGE', 'ISO', 'WYANT',
%                          'B&W', 'CIRCLE', 'HEXAGON', 'HEXAGONR30',
%                          'RECTANGLE', 'SQUARE', 'ANNULUS'
%
%   OUTPUT
%     results = The Zernike polynomial (n,m,sc) value for the given pixel (rho, theta).
%

% calculate radial part Rnm
Rnm = zeros(size(rho));
for s=0:(n-m)/2
    numerator = (-1)^s * factorial(n-s);
    denominator = factorial(s)*factorial((n+m)/2-s)*factorial((n-m)/2-s);
    Rnm = Rnm + (numerator / denominator) * rho.^(n-2*s);
end % for s statement

% 3 cases.  sc=' ', 's', 'c'.
theFactor = 1;
switch sc
    case ' '
        % means m=0
        theFactor = sqrt(n+1);
        result = Rnm;
    case 's'
        theFactor = sqrt(2*(n+1));
        result = Rnm .* sin(m*theta);
    case 'c'
        theFactor = sqrt(2*(n+1));
        result = Rnm .* cos(m*theta);
end % switch sc

switch ZernikeDef
    case {'MAHAJAN', 'NOLL', ...
            'HEXAGON', 'HEXAGONR30', 'RECTANGLE', 'SQUARE', 'ANNULUS'}
        result = theFactor * result;
end % switch


end % functon calculateZernike
