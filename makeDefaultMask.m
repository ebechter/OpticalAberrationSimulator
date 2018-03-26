function mask=makeDefaultMask(maskType, defaultRows, defaultCols, ShapeParameter)
% This function makes a default mask.  Since it is a default mask, the
% mask matrix is to be square.

mask = zeros(defaultRows, defaultCols);

% the circle into which the mask shape is to fit
cr = (defaultRows+1)/2;
cc = (defaultCols+1)/2;
defaultRadiusInPixels = (min(defaultRows, defaultCols) - 1)/2;

switch maskType
    case 'HEXAGON'
        % make a Hexagon mask
        for r=1:defaultRows
            for c=1:defaultCols
                x = (c-cc);
                y = (r-cr);
                rho = sqrt(x^2+y^2);
                eTheta = atan2(y,x);
                
                if (eTheta >= 0) && (eTheta <= (60*pi/180))
                    beta = 120*pi/180 - eTheta;
                    R = defaultRadiusInPixels * sind(60) / sin(beta);
                    if rho < R
                        mask(r,c) = 1;
                    end % if R statement
                end % if statement
                
                if eTheta >= (120*pi/180)
                    beta = 120*pi/180 - (pi - eTheta);
                    R = defaultRadiusInPixels * sind(60) / sin(beta);
                    if rho < R
                        mask(r,c) = 1;
                    end % if R statement
                end % if statement
                
                if eTheta <= -(120*pi/180)
                    beta = 120*pi/180 - (pi - abs(eTheta));
                    R = defaultRadiusInPixels * sind(60) / sin(beta);
                    if rho < R
                        mask(r,c) = 1;
                    end % if R statement
                end % if statement
                
                if (eTheta <= 0) && (eTheta >= -(60*pi/180))
                    beta = 120*pi/180 - abs(eTheta);
                    R = defaultRadiusInPixels * sind(60) / sin(beta);
                    if rho < R
                        mask(r,c) = 1;
                    end % if R statement
                end % if statement
                
                if (eTheta >= (60*pi/180)) && (eTheta <= (120*pi/180))
                    if abs(cr-r) <= ((sqrt(3)/2)*defaultRadiusInPixels)
                        mask(r,c) = 1;
                    end % if statement
                end % if statement
                
                if (eTheta <= (-60*pi/180)) && (eTheta >= (-120*pi/180))
                    if abs(cr-r) <= ((sqrt(3)/2)*defaultRadiusInPixels)
                        mask(r,c) = 1;
                    end % if statement
                end % if statement
                
                
            end % for c statement
        end % for r statement
        
        
    case 'HEXAGONR30'
        % make a Hexagon mask rotated 30 degrees
        for r=1:defaultRows
            for c=1:defaultCols
                x = (c-cc);
                y = (r-cr);
                rho = sqrt(x^2+y^2);
                eTheta = atan2(y,x);
                
                if (eTheta >= (30*pi/180)) && (eTheta <= (90*pi/180))
                    beta = 120*pi/180 - (eTheta - 30*pi/180);
                    R = defaultRadiusInPixels * sind(60) / sin(beta);
                    if rho < R
                        mask(r,c) = 1;
                    end % if R statement
                end % if statement
                
                if (eTheta >= (90*pi/180)) && (eTheta <= (150*pi/180))
                    beta = 120*pi/180 - (pi - eTheta - 30*pi/180);
                    R = defaultRadiusInPixels * sind(60) / sin(beta);
                    if rho < R
                        mask(r,c) = 1;
                    end % if R statement
                end % if statement
                
                if (eTheta <= -(90*pi/180)) && (eTheta >= -(150*pi/180))
                    beta = 120*pi/180 - (pi - abs(eTheta) - 30*pi/180);
                    R = defaultRadiusInPixels * sind(60) / sin(beta);
                    if rho < R
                        mask(r,c) = 1;
                    end % if R statement
                end % if statement
                
                if (eTheta <= -30*pi/180) && (eTheta >= -(90*pi/180))
                    beta = 120*pi/180 - (abs(eTheta) - 30*pi/180);
                    R = defaultRadiusInPixels * sind(60) / sin(beta);
                    if rho < R
                        mask(r,c) = 1;
                    end % if R statement
                end % if statement
                
                if (eTheta >= (-30*pi/180)) && (eTheta <= (30*pi/180))
                    if abs(cc-c) <= ((sqrt(3)/2)*defaultRadiusInPixels)
                        mask(r,c) = 1;
                    end % if statement
                end % if statement
                
                if (eTheta >= (150*pi/180)) || (eTheta <= (-150*pi/180))
                    if abs(cc-c) <= ((sqrt(3)/2)*defaultRadiusInPixels)
                        mask(r,c) = 1;
                    end % if statement
                end % if statement
                
                
            end % for c statement
        end % for r statement
        
    case 'ELLIPSE'
        % an ellipse mask is needed
        a = defaultRadiusInPixels;
        b = ShapeParameter * defaultRadiusInPixels;
        for r=1:defaultRows
            for c=1:defaultCols
                rho = sqrt((cr-r)^2 + (cc-c)^2);
                eTheta = atan2((cr-r),(cc-c));
                er = a*b/sqrt((b*cos(eTheta))^2+(a*sin(eTheta))^2);
                if rho <= er
                    mask(r,c) = 1;
                end % if statement
            end % for c statement
        end % for r stateemnt
    case 'RECTANGLE'
        % a rectangular mask is needed
        halfEdgec = ShapeParameter * defaultRadiusInPixels;
        halfEdger = sqrt(1 - ShapeParameter^2) * defaultRadiusInPixels;
        for r=1:defaultRows
            for c=1:defaultCols
                if (r > abs(cr-halfEdger)) && (r < (cr+halfEdger)) && ...
                        (c > abs(cr-halfEdgec)) && (c < (cr+halfEdgec))
                    mask(r,c) = 1;
                end % if statement
            end % for c statement
        end % for r stateemnt
    case 'SQUARE'
        % a square mask is needed
        halfEdge = (1/sqrt(2)) * defaultRadiusInPixels;
        for r=1:defaultRows
            for c=1:defaultCols
                if (r > abs(cr-halfEdge)) && (r < (cr+halfEdge)) && ...
                        (c > abs(cr-halfEdge)) && (c < (cr+halfEdge))
                    mask(r,c) = 1;
                end % if statement
            end % for c statement
        end % for r stateemnt
    case 'ANNULUS'
        % an annulus mask is needed
        obscuration = ShapeParameter * defaultRadiusInPixels;
        for r=1:defaultRows
            for c=1:defaultCols
                rho = sqrt((cr-r)^2 + (cc-c)^2);
                if (rho <= defaultRadiusInPixels) && (rho > obscuration)
                    mask(r,c) = 1;
                end % if statement
            end % for c statement
        end % for r stateemnt
    otherwise
        % a circle mask is needed
        for r=1:defaultRows
            for c=1:defaultCols
                rho = sqrt((cr-r)^2 + (cc-c)^2);
                if rho <= defaultRadiusInPixels
                    mask(r,c) = 1;
                end % if statement
            end % for c statement
        end % for r stateemnt
        
end % switch statement

end % function makeDefaultMask
