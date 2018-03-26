function ok = validateInput(ZernikeList, Zdata, mask, ...
    ZernikeDef, ShapeParameter, ...
    centerRow, centerCol, radiusInPixels, MMESorC)
%
%  Validates the input.
%  See above function definition for definition of input.
%
%   OUTPUT
%     ok  = true.
%

ok = true;

% Check validity of input parameters.

theIOCSize = size(Zdata);
theMaskSize = size(mask);

if (theIOCSize(1,2) > 1) && (sum(theIOCSize == theMaskSize) ~= 2)
    ME = MException('VerifyData:InvalidData', ...
        'ZernikeCalc: Surface data and mask matrices are not the same size.');
    throwAsCaller(ME);
end % if statement

if centerRow < 1
    ME = MException('VerifyData:InvalidData', ...
        'ZernikeCalc: centerRow must be positive.');
    throwAsCaller(ME);
end % if statement

if centerCol < 1
    ME = MException('VerifyData:InvalidData', ...
        'ZernikeCalc: centerCol must be positive.');
    throwAsCaller(ME);
end % if statement

if radiusInPixels < 1
    ME = MException('VerifyData:InvalidData', ...
        'ZernikeCalc: radiusInPixels must be positive.');
    throwAsCaller(ME);
end % if statement

hlda = mask == 0;
hldb = mask == 1;
if sum(sum(hlda + hldb)) ~= (theMaskSize(1,1)*theMaskSize(1,2))
    ME = MException('VerifyData:InvalidData', ...
        'ZernikeCalc: mask matrix must contain 0 or 1 only.');
    throwAsCaller(ME);
end % if statement


% Now for the fun stuff: Validating the Zernike ordering.

switch ZernikeDef
    case {'FRINGE', 'ISO', 'WYANT' 'MAHAJAN', 'NOLL', 'B&W', 'STANDARD', ...
            'HEXAGON', 'HEXAGONR30', 'ELLIPSE', 'RECTANGLE', ...
            'SQUARE', 'ANNULUS'}
        % These are the valid values.
    otherwise
        % ZernikeType is not valid
        ME = MException('VerifyData:InvalidData', ...
            'ZernikeCalc: ZernikeDef is not valid.');
        throwAsCaller(ME);
end % switch statement


theSize = size(ZernikeList);
rows = theSize(1,1);
cols = theSize(1,2);

hld1 = sum(sum(zeros(rows, cols) == (abs(ZernikeList) - floor(abs(ZernikeList)))));
if hld1 ~= rows*cols
    % a number in ZernikeList is not an integer
    ME = MException('VerifyData:InvalidData', ...
        'ZernikeCalc: ZernikeList must contain only integers.');
    throwAsCaller(ME);
end % if statement

if rows == 1
    % j ordering
    if sum(abs(ZernikeList(1, :)) ~= ZernikeList(1, :)) ~= 0
        ME = MException('VerifyData:InvalidData', ...
            'ZernikeCalc: ZernikeList j values must be positive or 0.');
        throwAsCaller(ME);
    end % if statement
    
    switch ZernikeDef
        case 'FRINGE'
            if sum(ZernikeList(1, :) < 1) ~= 0
                ME = MException('VerifyData:InvalidData', ...
                    'ZernikeCalc: FRINGE order value j must be greater than 0.');
                throwAsCaller(ME);
            end % if statement
            
            if sum(ZernikeList(1, :) > 37) ~= 0
                ME = MException('VerifyData:InvalidData', ...
                    'ZernikeCalc: FRINGE order value j must be less than 38.');
                throwAsCaller(ME);
            end % if statement
            
        case 'ISO'
            if sum(ZernikeList(1, :) > 35) ~= 0
                ME = MException('VerifyData:InvalidData', ...
                    'ZernikeCalc: ISO order value j out of range.');
                throwAsCaller(ME);
            end % if statement
            
        case 'WYANT'
            if sum(ZernikeList(1, :) > 36) ~= 0
                ME = MException('VerifyData:InvalidData', ...
                    'ZernikeCalc: WYANT order value j out of range.');
                throwAsCaller(ME);
            end % if statement
            
        case {'MAHAJAN', 'NOLL', 'B&W', 'STANDARD', ...
                'HEXAGON', 'HEXAGONR30', 'ELLIPSE', 'RECTANGLE', ...
                'SQUARE', 'ANNULUS'}
            if sum(ZernikeList(1, :) == 0) ~= 0
                ME = MException('VerifyData:InvalidData', ...
                    'ZernikeCalc: Zernike order value j can not be zero.');
                throwAsCaller(ME);
            end % if statement
            
    end % switch statement
    
else
    % (n, m) specification
    if sum(abs(ZernikeList(1, :)) ~= ZernikeList(1, :)) ~= 0
        ME = MException('VerifyData:InvalidData', ...
            'ZernikeCalc: ZernikeList n values must be positive or 0.');
        throwAsCaller(ME);
    end % if statement
    
    if sum(abs(ZernikeList(2, :)) > ZernikeList(1, :)) ~= 0
        ME = MException('VerifyData:InvalidData', ...
            'ZernikeCalc: abs(m) values must not exceed n values.');
        throwAsCaller(ME);
    end % if statement
    
    switch ZernikeDef
        case 'FRINGE'
            if rows > 1
                ME = MException('VerifyData:InvalidData', ...
                    'ZernikeCalc: FRINGE can only be specified with j ordering.');
                throwAsCaller(ME);
            end % if statement
            
        case 'ISO'
            if rows > 1
                ME = MException('VerifyData:InvalidData', ...
                    'ZernikeCalc: ISO can only be specified with j ordering.');
                throwAsCaller(ME);
            end % if statement
            
        case 'WYANT'
            if rows > 1
                ME = MException('VerifyData:InvalidData', ...
                    'ZernikeCalc: WYANT can only be specified with j ordering.');
                throwAsCaller(ME);
            end % if statement
            
        case {'HEXAGON', 'HEXAGONR30', 'ELLIPSE', 'RECTANGLE', ...
                'SQUARE', 'ANNULUS'}
            if rows > 1
                ME = MException('VerifyData:InvalidData', ...
                    'ZernikeCalc: ZernikeDef value can only be specified with j ordering.');
                throwAsCaller(ME);
            end % if statement
            
        otherwise
            if (sum(mod(ZernikeList(1, :) - abs(ZernikeList(2, :)), 2)) ~= 0)
                ME = MException('VerifyData:InvalidData', ...
                    'ZernikeCalc: n-abs(m) must be even.');
                throwAsCaller(ME);
            end % if statement
    end % switch
    
    switch MMESorC
        case {'-m=sin', '-m=cos'}
        otherwise
            ME = MException('VerifyData:InvalidData', ...
                'ZernikeCalc: Sign convention value for (n, m) not valid.');
            throwAsCaller(ME);
    end % switch
    
end % (n,m) order specified

% Validate ZernikeShape parameter

switch ZernikeDef
    case {'ELLIPSE', 'RECTANGLE','ANNULUS'}
        % Check that the ShapeParameter is valid
        if (ShapeParameter <= 0) || (ShapeParameter >= 1)
            ME = MException('VerifyData:InvalidData', ...
                'ZernikeCalc: The ShapeParameter is out of range: (0,1).');
            throwAsCaller(ME);
        end % if statement
end % switch ZernikeDef

end % function validateInput
