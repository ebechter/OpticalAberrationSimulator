function [output1 output2] = ZernikeCalc( ZernikeList, Zdata, mask, ...
    ZernikeDef, ShapeParameter, ...
    unitCircle, MMESorC)
%ZERNIKECALC Uses 'mask' region to fit circular (or other shape) Zernikes to surface data.
%
% VERSION:  2013-02-07   (YYYY-MM-DD)
%
% Fits circular, hexagonal, rectangular, square, elliptical, or annulus
% orthogonal Zernike polynomials to the surface data provided.  If no surface
% data is provided then plots of these polynomial functions over the
% region specified are produced.  Several different convesions (sign,
% normalization) can be explicitly specified.
%
% function [output1 output2] = ZernikeCalc( ZernikeList, Zdata, mask, ...
%                                           ZernikeDef, ShapeParameter, ...
%                                           unitCircle, MMESorC)
%
%  INPUT:
%    ZernikeList     = if size() = (1,k) then this is a Zernike j-ordering
%                      list of numbers (j).
%                      if size() = (2,k) then this is Zernike (n, m) list.
%                      row 1 is the list of n values, row 2 is the list of
%                      m values.  m is positive or negative values to indicate
%                      whether sin() or cos() factor for the term.   See MMESorC
%                      parameter for the meaning of a minus m value.
%                      DEFAULT: [1:15]. First 15 j-ordered Zernikes.
%
%    Zdata           = If this is a single column then it contains the
%                      Zernike coefficients and we are to calculate the
%                      surface data.
%                      If this is more than a single column then it is a
%                      surface data and we are to calculate the
%                      coefficients for the Zernike polynomials specified.
%                      The surface data must be a phase unwrapped surface data.
%                      DEFAULT: ones(15,1).
%
%    mask            = Matrix (same size as surface data) indicating
%                      the individual pixels in the surface data that are
%                      to be used for fitting of the Zernike polynomials.
%                      '0' = don't use this pixel, '1' = use this pixel.
%                      If mask is a scalar number n, then an nxn mask matrix
%                      is created with an n diameter circle as the mask.
%                      DEFAULT: 100x100 with circluar mask.
%
%    ZernikeDef      = One of 'FRINGE', 'ISO', 'WYANT', 'MAHAJAN', 'NOLL',
%                      'B&W', 'STANDARD'
%                      'HEXAGON', 'HEXAGONR30', 'ELLIPSE',
%                      'RECTANGLE', 'SQUARE', 'ANNULUS'
%                      See table below for possible j-value ranges.
%                      NOTE: 'MAHAJAN' and 'NOLL' are programmed to be the same.
%                      NOTE: 'B&W' and 'STANDARD' are programmed to be the same.
%                      DEFAULT: 'FRINGE'.
%
%    ShapeParameter  = For 'ELLIPSE' this is the aspect ratio (b).
%                      For 'RECTANGLE' this is half-width along the x-axis (a).
%                      For 'ANNULUS' this is the obscuration ratio (e).
%                      For all other shapes this is ignored.
%                      DEFAULT: 0, which is not valid for 'ELLIPSE' and 'RECTANGLE'.
%
%    unitCircle      = Gives the unit circle around the mask which the circular
%                      Zernikes are defined on.  It is at most a 1x3 vector
%                      specifying centerRow, centerCol, radiusInPixels.  The
%                      permutations and ording allowed are as follows:
%                      []  empty so all 3 defaults are calculated.
%                      [radiusInPixels]  the default centerRow and centerCol
%                                        are calculated.
%                      [centerRow, centerCol] the default radiusInPixels is
%                                        calculated.
%                      [centerRow, centerCol, radiusInPixels]  no defaults.
%                      DEFAULT: The defaults for centerRow and centerCol are
%                               calculated by  calculating the weighted
%                               average of row and column for the mask's 1's.
%                               The default radiusInPixels is calculated to
%                               be the largest distance from (centerRow,centerCol)
%                               to all mask's pixels with a '1' value.
%
%    MMESorC         = Either '-m=sin' or '-m=cos'.  Indicates, for (n, m) ordering
%                      what a minus m value represents, a sine factor or a
%                      cosine factor.
%                      DEFAULT: '-m=sin'. Ignored when j ordering is used.
%
%  OUTPUT:
%    - When no output result is requested (ZernikeCalc(...)) this function
%      produces plots of the Zernike polynomials requested.
%      If Zdata is a single column specifying the coefficients then
%          a plot of the sum of each of the Zernikes specified is produced.
%      If Zdata is an n x m matrix of surface data then 3 plots are
%      produced: 1) of the surface data, 2) of the fitted surface, 3) the
%      difference between the surface data and the fitted surface.
%    - When one output result is requested (results = ZernikeCalc(...)) then 1
%      of 2 possible results are returned:
%        if the input Zdata is a single column specifying the coefficients
%           to multiply the Zernike polynomials by then the result retruned
%           is the Zernike polynomial value matrices (across the mask area).
%        if the input Zdata is a matrix for which the Zernikes are to be
%           fit then the results returned is a column vector of the Zernike
%           coefficients corrseponding to (and in the same order as) the
%           polynomials identified by ZernikeList.
%    - If 2 output results are requested ([out1 out2] = ZernikeCalc(...)) then
%            if the input for Zdata is a single column, giving the coefficients
%               to multiply the Zernike polynomials by, then out1 is the sum
%               of the Zernike polynomials requested, and out2 is a 3 dimensional
%               matrix of all the n Zernike polynomials requested [:,:,n].
%            if the input for Zdata is not a column vector then out1 is the
%                3 dimensional data cube of the n Zernike polynomials requested [:,:,n]
%                and out2 are the coefficients used.
%
%
% Examples:
%
%   ZernikeCalc
%   - Displays the first 15 Fringe Zernikes in 15 color plots.
%
%   ZernikeCalc([4 5 7 8], [0.4; -0.6; 1.2; 0.25])
%   - Displays Fringe Zernikes Zj=4, Zj=5, Zj=7, Zj=8 multiplied by the
%   scaling factors 0.4, -0.6, 1.2 and 0.25, respectively in 4 separate
%   color plots.
%
%   ZernikeCalc([4 5 7 8], [0.4; -0.6; 1.2; 0.25], [], 'standard')
%   - Same as the last case only using standard Zernikes rather than Fringe
%   Zernikes.
%
%   ZernikeCalc([2 2; 2 0; 3 3; 3 1]', [0.4; -0.6; 1.2; 0.25], [], 'standard')
%   - Same as last case now using Z(n,m) notation to specify which Zernikes
%   to use.
%
%   Let SD be an n x m matrix of surface data to which the specified
%   Zernikes are to be fit.  Then
%
%   coeff = ZernikeCalc([2 2; 2 0; 3 3; 3 1]', SD, [], 'standard')
%
%   returns a column consisting of the calculated fitting coefficients
%   in coeff.  No plots are produced.
%
%   [DC, coeff] = ZernikeCalc([2 2; 2 0; 3 3; 3 1]', SD, [], 'standard')
%
%   returns a column consisting of the calculated fitting coefficients
%   in coeff and a n x m x 4 data cube.  ('4' because 4 Zernikes were
%   specified) Each DC(:, :, i) (i=1,2,3,4) is the ith specified Zernike
%   fitted to the surface data SD across the (default) mask area.
%
%   [DC, coeff] = ZernikeCalc([4 5 7 8], SD, [], 'annulus', 0.25)
%
%   This uses the annular Zernikes, with a central obscuration radius ratio of 0.25
%   to fit the surface data.  See Ref. 1 for details on noncircular Zernikes.
%
%
% Circular Zernike polynomials are available in several optical design
% software packages, including Code V(R), OSLO(R), Zemax(R), etc.
%
%      Table 1 Software Conventions
%--------------------------------------------
%|INPUT PARAM |APERTURE |SOFTWARE|ORDER &   |
%|ZernikeDef  | SHAPE   |        | RANGE    |
%|------------|---------|--------|----------|
%|'B&W',      |Circle   |CODE V  | (n, ±m), |
%|'STANDARD'  |         |        | j=1...   |
%|------------|---------|--------|----------|
%|'MAHAJAN',  |Circle   |ZEMAX   | (n, ±m), |
%|'NOLL'      |         |        | j=1...   |
%|------------|---------|--------|----------|
%|'FRINGE'    |Circle   |CODE V, | j=1...37 |
%|            |         |ZEMAX   |          |
%|------------|---------|--------|----------|
%|'ISO'       |Circle   |        | j=0..35  |
%|------------|---------|--------|----------|
%|'WYANT'     |Circle   | OSLO   | j=0...36 |
%|------------|---------|--------|----------|
%|'HEXAGON'   |Hexagon  |        | j=1...45 |
%|------------|---------|--------|----------|
%|'HEXAGONR30'|Heaxgon  |        | j=1...28 |
%|            |rotated  |        |          |
%|            |30 deg.  |        |          |
%|------------|---------|--------|----------|
%|'ELLIPSE'   |Ellipse  |        | j=1..15  |
%|------------|---------|--------|----------|
%|'RECTANGLE' |Rectangle|        | j=1...15 |
%|------------|---------|--------|----------|
%|'SQUARE'    |Square   |        | j=1...45 |
%|------------|---------|--------|----------|
%|'ANNULUS'   |Annulus  |        | j=1...35,|
%|            |         |        | j<>29,30,|
%|            |         |        |    31,32 |
%--------------------------------------------
%
% Ref. 1:  Mahajan, V.N., G.-m. Dai, "Orthonormal polynomials in wavefront
%          analysis: analytical solution," J. Opt. Soc. Am. A, Vol. 24, No. 9
%          Sept. 2007.
%

%
% Updates:  2012-01-08  (YYYY-MM-DD)
%           RWG - Added default mask shapes for the different ZernikeDef
%                 input parameter values.
%
%           2012-01-08  (YYYY-MM-DD)
%           RWG - When no output requested ZernikeCalc will print all
%                 Zernike polynomials specified.
%


%
% Code Copyrighted, 2011-2013 by Robert W. Gray and Joseph M. Howard.  All
% rights reserved.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For validation of the equations, uncomment the next 2 lines.
% Then, it doesn't matter what input parameters are specified.
%

%          validateZ();
%          return;

%
% end validation of equations.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assign default values to input parameters that are empty.
% Empty is '' for strings and [] for vectors and matrices.
%

if nargin == 0 || isempty(ZernikeList)
    ZernikeList = 1:15;    % default is first 15 fringe Zernikes.
end % if statement

if nargin <= 1 || isempty(Zdata)
    theSize = size(ZernikeList);
    Zdata = ones(theSize(1,2),1);  % all coefficients are 1.
end % if statement


if nargin <= 3 || isempty(ZernikeDef)
    ZernikeDef = 'FRINGE';
end % if statement

% Convert to upper case
ZernikeDef = upper(ZernikeDef);


if nargin <= 4 || isempty(ShapeParameter)
    ShapeParameter = 0;
end % if statement


if nargin <= 2 || isempty(mask)
    % make a default mask
    
    defaultRows = 100;
    defaultCols = 100;
    
    theZdataSize = size(Zdata);
    if (theZdataSize(1,1) > 1) && (theZdataSize(1,2) > 1)
        defaultRows = theZdataSize(1,1);
        defaultCols = theZdataSize(1,2);
    end % if statement
    
    mask = makeDefaultMask(ZernikeDef, defaultRows, defaultCols, ShapeParameter);
    
end % if no mask statement

sm = size(mask);
if (sm(1,1) == 1) && (sm(1,2) == 1)
    % if mask is a scalar n then create nxn matrix
    % make a circular mask
    defaultRows = mask(1,1);
    defaultCols = mask(1,1);
    
    mask = makeDefaultMask(ZernikeDef, defaultRows, defaultCols, ShapeParameter);
    
end % end if statement




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate default centerRow and centerCol values.
% These might have been specified, but if they haven't
% use these values.

radiusInPixels = 0;

theSize = size(mask);
numrows = theSize(1,1);
numcols = theSize(1,2);

% calculate center by weighted averages.
sumMaskRow = sum(mask,2);
sumMaskCol = sum(mask,1);

sumMaskAll = sum(sum(mask));

centerRow = sum(sumMaskRow .* ((1:numrows)')) / sumMaskAll;
centerCol = sum(sumMaskCol .* (1:numcols)) / sumMaskAll;

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin <= 5 || isempty(unitCircle)
    % The unit circle associated with the mask has not been specified.
    % unitCircle is empty.
    
    unitCircle = [centerRow, centerCol];
    
end % if statement

% At this point, unitCircle is not empty.
theSize = size(unitCircle);
nc = theSize(1,2);
switch nc
    case 1
        % only the radiusInPixels has been specified
        radiusInPixels = unitCircle(1,1);
        
    case 2
        % the centerRow and centerCol have been specified
        % so can now calculate the radius in Pixels.
        centerRow = unitCircle(1,1);
        centerCol = unitCircle(1,2);
        
        % a matrix such that each element (r,c) has the value r-centerRow
        rm = (((1:numrows)-centerRow)'*ones(1,numcols)).*mask;
        % a matrix such that each element (r,c) has the value c-centerCol
        cm = (ones(numrows,1)*((1:numcols)-centerCol)).*mask;
        % sqrt(rm.^2 + cm.^2) is a matrix such that (r,c) contains the distance
        % of (r,c) to the center (centerRow, centerCol).
        radiusInPixels = max(max(sqrt(cm.^2 + rm.^2)));
        
    case 3
        % the centerRow, centerCol, radiusInPixels have been specified
        centerRow = unitCircle(1,1);
        centerCol = unitCircle(1,2);
        radiusInPixels = unitCircle(1,3);
        
    otherwise
        % error.
end % switch statement

if nargin <= 6 || isempty(MMESorC)
    MMESorC = '-m=sin';
end % if stateemnt

%
% end of section on default input assignments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input validation
% Too much code is distracting, so put validation in separate function.
validateInput(ZernikeList, Zdata, mask,  ...
    ZernikeDef, ShapeParameter, ...
    centerRow, centerCol, radiusInPixels,  MMESorC);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The Zernike polynomials are caluclated using (n, m) ordering
% so need to calculate (n,m,sc) even when j ordering is specified.
%  sc = 's' means sin() factor
%  sc = 'c' means cos() factor
%  sc = ' ' means m = 0 so has no sin() or cos() factor
%

theSize = size(ZernikeList);
maxNumOfZernikeTerms = theSize(1,2);
n = zeros(1, theSize(1,2));
m = zeros(1, theSize(1,2));

sc = ' ';
if theSize(1,1) == 1
    % the ZernikeList is list of j order values.
    for k=1:theSize(1,2)
        % convert the j ordering of type ZernikeDef to (n,m,sc) of
        % same ZernikeDef.
        [n(k) m(k) sc(k)] = convertj(ZernikeList(1, k), ZernikeDef);
    end % for statement
else
    % the ZernikeList is list of (n,m) pairs.
    for k=1:theSize(1,2)
        % convert (n,m) to (n,m,sc) using MMESorC
        n(k) = ZernikeList(1, k);
        m(k) = abs(ZernikeList(2, k));
        sc(k) = ' ';
        switch MMESorC
            case '-m=sin'
                if ZernikeList(2, k) < 0
                    sc(k) = 's';
                end % if statement
                if ZernikeList(2, k) > 0
                    sc(k) = 'c';
                end % if statement
            case '-m=cos'
                if ZernikeList(2, k) < 0
                    sc(k) = 'c';
                end % if statement
                if ZernikeList(2, k) > 0
                    sc(k) = 's';
                end % if statement
        end % switch statement
        
    end % for k statement
    
end % if j or (n,m) ordering

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preallocate some of the matrices.
%

% We'll need the size of the mask matrix (same as surface data).
theSize = size(mask);
rows = theSize(1,1);
cols = theSize(1,2);

% pre-allocate vectors and matrices.
numOfMaskPixels = sum(sum(mask));

zcoeff = zeros(maxNumOfZernikeTerms, 1);
xMatrices = zeros(rows, cols, maxNumOfZernikeTerms);
xVectors = zeros(numOfMaskPixels, maxNumOfZernikeTerms);
yVector = zeros(numOfMaskPixels, 1);

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check if surface data is the input or if a coefficent vector is
% the input.

interferogramInput = true;

theSize = size(Zdata);
if theSize(1,2) == 1
    % There is one column in Zdata so it is a column vector of coefficents.
    zcoeff = Zdata;
    interferogramInput = false;
end % if statement

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  for each pixel (matrix element) we need the distance from the
%  center (centerRow, centerCol).  So we make a matrix that has
%  as element values the distance of that (row, col) element to
%  (centerRow, centerCol).
%

% a matrix such that each element (r,c) has the value r-centerRow
rm = ((1:rows)-centerRow)'*ones(1,cols);
% a matrix such that each element (r,c) has the value c-centerCol
cm = ones(rows,1)*((1:cols)-centerCol);

% sqrt(rm.^2 + cm.^2)./radiusInPixels is a matrix such that (r,c) contains
% the normalized distance of (r,c) to the center (centerRow, centerCol).
% Then we reshape this into a column vector.
rho = reshape(sqrt(cm.^2 + rm.^2)./radiusInPixels, rows*cols, 1);
% atan2(rm, cm) is a matrix such that each element (r,c) contains
% the angle (radians) from the centerRow axis to (r,c).  Then we
% reshape this into a vector.
theta = reshape(atan2(rm, cm), rows*cols, 1);
% reshape the mask into a vector.
vmask = reshape(mask, rows*cols, 1);

if interferogramInput
    % we have surface data so reshape it just like rho and theta.
    yVector = reshape(Zdata, rows*cols, 1);
    % and remove all the elements that don't have a corresponding '1'
    % in the mask.
    yVector(vmask ==0) = [];
end % if statement

for i=1:maxNumOfZernikeTerms
    % calculate the ith Zernike polynomial value for each pixel
    % in the mask matrix (now a vector).
    
    switch ZernikeDef
        % Handle the special shapes.
        case 'HEXAGON'
            hldZ = ZHexagon(ZernikeList(1,i), rho, theta);
        case 'HEXAGONR30'
            hldZ = ZHexagonR30(ZernikeList(1,i), rho, theta);
        case 'ELLIPSE'
            hldZ = ZEllipse(ZernikeList(1,i), rho, theta, ShapeParameter);
        case 'RECTANGLE'
            hldZ = ZRectangle(ZernikeList(1,i), rho, theta, ShapeParameter);
        case 'SQUARE'
            hldZ = ZSquare(ZernikeList(1,i), rho, theta);
        case 'ANNULUS'
            hldZ = ZAnnulus(ZernikeList(i), n(i), m(i), rho, theta, ShapeParameter);
        otherwise
            % Otherwise, its a circle Zernike.
            hldZ = calculateZernike(n(i), m(i), sc(i), rho, theta, ZernikeDef);
            
    end % switch ZernikeShape
    
    % reshape the column vector result into a (rows, cols) matrix.
    xMatrices(:, :, i) = reshape(hldZ,rows,cols);
    
    % remove the elements from the Zernike calculation that do not
    % have a corresponding '1' in the mask.
    hldZ(vmask == 0) = [];
    % this is one of the Zernike polynomial results for each pixel
    % for which the mask is a '1'.
    xVectors(:, i) = hldZ;
    
end % for i statement

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do the regression (least squares fit) only if the
% input Zdata is surface data.
%

if interferogramInput
    % Use least squares fit to determine the coefficients.
    zcoeff = xVectors\yVector;
end % if statement
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% multiply Zernike polynomial matrices by the calculated coefficients.
% since we already have the Zernike matrices calculated, we do this here.
%
zm = zeros(rows, cols, maxNumOfZernikeTerms);
for i=1:maxNumOfZernikeTerms
    % multiply the Zernike matrix by the corresponding coefficient
    % and by the mask to zero out pixels that we don't care about.
    zm(:,:,i) =  xMatrices(:, :, i) .* zcoeff(i) .* mask;
end % for statement
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Change the output depending on input and output specified.
%
if nargout < 1
    % plot Zernike figures
    sizeZD = size(Zdata);
    
    if sizeZD(1,2) == 1
        % plot each of the Zernike polynomials
        figs(zm);
        % plot the sum of the request Zernike Polynomials
        %       figs(sum(zm,3));
        return;
    end % if statement
    
    % plot 1) input surface data, 2) the fit results, 3) their difference
    
    toPlot(:,:,1) = Zdata;
    
    theFit = sum(zm, 3);
    toPlot(:,:,2) = theFit;
    
    theDiff = Zdata - theFit;
    toPlot(:,:,3) = theDiff;
    
    figs(toPlot, mask, 3, {'Input Data', 'Zernike Fit', 'Input Data - Zernike Fit'});
    
    return;
end % if nargout < 1

% At this point we know nargout >= 1

if nargout < 2
    % This means nargout == 1.
    if interferogramInput
        % return calculated coefficients only
        output1 = zcoeff;
        return;
    end % if interferogramInput
    
    % Not an interferogram as Input.  Must be coefficients as input
    % so return only the interferogram.  They already know the
    % coefficients.
    
    output1 = zm;
    return;
end % if nargout < 2

% At this point we know the output requested is 2 (or more).

theZdataSize = size(Zdata);
if theZdataSize(1,2) == 1
    output1 = sum(zm,3);
    output2 = zm;
    return;
end % if statement

output1 = zm;
output2 = zcoeff;

end % ZernikeCalc