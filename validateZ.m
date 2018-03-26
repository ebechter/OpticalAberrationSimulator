function validateZ()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%  Validate Circle

%     fprintf(1,'\n\n    CIRCLE  VALIDATION \n');
%
%     ZernikeType = 'MAHAJAN';
%
%     for j1=1:37
%
%         j2 = j1 - 1;
%         if j1 == 1
%           j2 = 2;
%         end % if statement
%
%         [n1,m1,sc1] = convertj(j1, ZernikeType);
%         [n2,m2,sc2] = convertj(j2, ZernikeType);
%
%         fun1 = @(rho,theta) calculateZernike(n1,m1,sc1,rho,theta,ZernikeType).*calculateZernike(n1,m1,sc1,rho,theta,ZernikeType).*rho;
%         Q1 = quad2d(fun1, 0,1, 0,2*pi);
%         diff1 = pi - Q1;
%
%         fun2 = @(rho,theta) calculateZernike(n1,m1,sc1,rho,theta,ZernikeType).*calculateZernike(n2,m2,sc2,rho,theta,ZernikeType).*rho;
%         Q2 = quad2d(fun2, 0,1, 0,2*pi);
%         diff2 = 0.0 - Q2;
%
%         fprintf(1,'j1=%d, n1=%d, m1=%d, result = %+20.18f  Difference = %+20.18f \n',j1,n1,m1,Q1,diff1);
%         fprintf(1,'j2=%d, n2=%d, m2=%d, result = %+20.18f  Difference = %+20.18f \n\n',j2,n2,m2,Q2,diff2);
%
%
%     end % for statement

%
% end Circle validation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%  Validate Circle Fringe

%     fprintf(1,'\n\n    CIRCLE (FRINGE)  VALIDATION \n');
%
%     ZernikeType = 'FRINGE';
%
%     for j1=1:37
%
%         j2 = j1 - 1;
%         if j1 == 1
%           j2 = 2;
%         end % if statement
%
%         [n1,m1,sc1] = convertj(j1, ZernikeType);
%         [n2,m2,sc2] = convertj(j2, ZernikeType);
%
%         fun1 = @(rho,theta) calculateZernike(n1,m1,sc1,rho,theta,ZernikeType).*calculateZernike(n1,m1,sc1,rho,theta,ZernikeType).*rho;
%         Q1 = quad2d(fun1, 0,1, 0,2*pi);
%         % The normalization is different.
%         em = 1;
%         if m1 == 0
%           em = 2;
%         end % if statement
%         diff1 = (em*pi)/(2*n1+2) - Q1;
%
%         fun2 = @(rho,theta) calculateZernike(n1,m1,sc1,rho,theta,ZernikeType).*calculateZernike(n2,m2,sc2,rho,theta,ZernikeType).*rho;
%         Q2 = quad2d(fun2, 0,1, 0,2*pi);
%         diff2 = 0.0 - Q2;
%
%         fprintf(1,'j1=%d, n1=%d, m1=%d, result = %+20.18f  Difference = %+20.18f \n',j1,n1,m1,Q1,diff1);
%         fprintf(1,'j2=%d, n2=%d, m2=%d, result = %+20.18f  Difference = %+20.18f \n\n',j2,n2,m2,Q2,diff2);
%
%
%     end % for statement

%
% end Circle Fringe validation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%  Validate ZHexagon

%     fprintf(1,'\n\n    HEXAGON  VALIDATION \n');
%     for j1=1:45
%
%         % Need another j value to compare orthogonality condition
%         j2 = j1 - 1;
%         if j1 == 1
%           j2 = 2;
%         end % if statement
%
%
%         ymina = @(x) -(sqrt(3)*x + sqrt(3));
%         ymaxa = @(x)   sqrt(3)*x + sqrt(3);
%
%         yminc = @(x) -(-sqrt(3)*x + sqrt(3));
%         ymaxc = @(x)   -sqrt(3)*x + sqrt(3);
%
%         fun1 = @(x,y) ZHexagon(j1,sqrt(x.*x+y.*y),atan2(y,x)).*ZHexagon(j1,sqrt(x.*x+y.*y),atan2(y,x));
%         Q1a = quad2d(fun1, -1,-0.5, ymina,ymaxa);
%         Q1b = quad2d(fun1, -0.5,0.5, -sqrt(3)/2,sqrt(3)/2);
%         Q1c = quad2d(fun1, 0.5,1, yminc,ymaxc);
%         diff1 = 3*sqrt(3)/2 - (Q1a+Q1b+Q1c);
%
%         fun2 = @(x,y) ZHexagon(j1,sqrt(x.*x+y.*y),atan2(y,x)).*ZHexagon(j2,sqrt(x.*x+y.*y),atan2(y,x));
%         Q2a = quad2d(fun2, -1,-0.5, ymina,ymaxa);
%         Q2b = quad2d(fun2, -0.5,0.5, -sqrt(3)/2,sqrt(3)/2);
%         Q2c = quad2d(fun2, 0.5,1, yminc,ymaxc);
%         diff2 = 0.0 - (Q2a+Q2b+Q2c);
%
%         fprintf(1,'j1=%d, result = %+20.18f  Difference = %+20.18f \n',j1,(Q1a+Q1b),diff1);
%         fprintf(1,'j2=%d, result = %+20.18f  Difference = %+20.18f \n\n',j2,(Q2a+Q2b),diff2);
%
%
%     end % for statement

%
% end ZHexagon validation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%  Validate ZHexagonR30

%     fprintf(1,'\n\n    HEXAGON ROTATED 30 DEGREES  VALIDATION \n');
%     for j1=1:28
%
%         % Need another j value to compare orthogonality condition
%         j2 = j1 - 1;
%         if j1 == 1
%           j2 = 2;
%         end % if statement
%
%
%         ymina = @(x) -(x/sqrt(3) + 1);
%         ymaxa = @(x)   x/sqrt(3) + 1;
%
%         yminb = @(x) -(-x/sqrt(3) + 1);
%         ymaxb = @(x)   -x/sqrt(3) + 1;
%
%         fun1 = @(x,y) ZHexagonR30(j1,sqrt(x.*x+y.*y),atan2(y,x)).*ZHexagonR30(j1,sqrt(x.*x+y.*y),atan2(y,x));
%         Q1a = quad2d(fun1, -sqrt(3)/2,0, ymina,ymaxa);
%         Q1b = quad2d(fun1, 0,sqrt(3)/2, yminb,ymaxb);
%         diff1 = 3*sqrt(3)/2 - (Q1a+Q1b);
%
%         fun2 = @(x,y) ZHexagonR30(j1,sqrt(x.*x+y.*y),atan2(y,x)).*ZHexagonR30(j2,sqrt(x.*x+y.*y),atan2(y,x));
%         Q2a = quad2d(fun2, -sqrt(3)/2,0, ymina,ymaxa);
%         Q2b = quad2d(fun2, 0,sqrt(3)/2, yminb,ymaxb);
%         diff2 = 0.0 - (Q2a+Q2b);
%
%         fprintf(1,'j1=%d, result = %+20.18f  Difference = %+20.18f \n',j1,(Q1a+Q1b),diff1);
%         fprintf(1,'j2=%d, result = %+20.18f  Difference = %+20.18f \n\n',j2,(Q2a+Q2b),diff2);
%
%
%     end % for statement

%
% end ZHexagonR30 validation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%  Validate ZEllipse
%
%     b = 0.8;  % semi-minor axis
%     %b = 1.0;  % semi-minor axis
%
%     fprintf(1,'\n\n    ELLIPSE  VALIDATION (b=%f)\n',b);
%     for j1=1:15
%
%         % Need another j value to compare orthogonality condition
%         j2 = j1 - 1;
%         if j1 == 1
%           j2 = 2;
%         end % if statement
%
%
%         ymax = @(x) b*sqrt(1-x.*x);
%         ymin = @(x) -b*sqrt(1-x.*x);
%
%         fun1 = @(x,y) ZEllipse(j1,sqrt(x.*x+y.*y),atan2(y,x),b).*ZEllipse(j1,sqrt(x.*x+y.*y),atan2(y,x),b);
%         Q1 = quad2d(fun1, -1,1, ymin,ymax);
%         diff1 = pi*b - Q1;
%
%         fun2 = @(x,y) ZEllipse(j1,sqrt(x.*x+y.*y),atan2(y,x),b).*ZEllipse(j2,sqrt(x.*x+y.*y),atan2(y,x),b);
%         Q2 = quad2d(fun2, -1,1, ymin,ymax);
%         diff2 = 0.0 - Q2;
%
%         fprintf(1,'j1=%d, result = %+20.18f  Difference = %+20.18f \n',j1,Q1,diff1);
%         fprintf(1,'j2=%d, result = %+20.18f  Difference = %+20.18f \n\n',j2,Q2,diff2);
%
%
%     end % for statement

%
% end ZEllipse validation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%  Validate ZRectangle

%     a = 0.37;  % half length of rectangle
%     %a = 1/sqrt(2);  % for a square
%
%     fprintf(1,'\n\n    RECTANGLE  VALIDATION (a=%f)\n',a);
%     for j1=1:15
%
%         % Need another j value to compare orthogonality condition
%         j2 = j1 - 1;
%         if j1 == 1
%           j2 = 2;
%         end % if statement
%
%         fun1 = @(x,y) ZRectangle(j1,sqrt(x.*x+y.*y),atan2(y,x),a).*ZRectangle(j1,sqrt(x.*x+y.*y),atan2(y,x),a);
%         Q1 = quad2d(fun1, -a,a, -sqrt(1-a*a),sqrt(1-a*a));
%         diff1 = (2*a * 2*sqrt(1-a*a)) - Q1;
%
%         fun2 = @(x,y) ZRectangle(j1,sqrt(x.*x+y.*y),atan2(y,x),a).*ZRectangle(j2,sqrt(x.*x+y.*y),atan2(y,x),a);
%         Q2 = quad2d(fun2, -a,a, -sqrt(1-a*a),sqrt(1-a*a));
%         diff2 = 0.0 - Q2;
%
%         fprintf(1,'j1=%d, result = %+20.18f  Difference = %+20.18f \n',j1,Q1,diff1);
%         fprintf(1,'j2=%d, result = %+20.18f  Difference = %+20.18f \n\n',j2,Q2,diff2);
%
%
%     end % for statement

%
% end ZRectangle validation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%  Validate ZSquare

%     fprintf(1,'\n\n    SQUARE  VALIDATION \n');
%     for j1=1:45
%
%         % Need another j value to compare orthogonality condition
%         j2 = j1 - 1;
%         if j1 == 1
%           j2 = 2;
%         end % if statement
%
%         fun1 = @(x,y) ZSquare(j1,sqrt(x.*x+y.*y),atan2(y,x)).*ZSquare(j1,sqrt(x.*x+y.*y),atan2(y,x));
%         Q1 = quad2d(fun1, -1/sqrt(2),1/sqrt(2), -1/sqrt(2),1/sqrt(2));
%         diff1 = (sqrt(2)*sqrt(2)) - Q1;
%
%         fun2 = @(x,y) ZSquare(j1,sqrt(x.*x+y.*y),atan2(y,x)).*ZSquare(j2,sqrt(x.*x+y.*y),atan2(y,x));
%         Q2 = quad2d(fun2, -1/sqrt(2),1/sqrt(2), -1/sqrt(2),1/sqrt(2));
%         diff2 = 0.0 - Q2;
%
%         fprintf(1,'j1=%d, result = %+20.18f  Difference = %+20.18f \n',j1,Q1,diff1);
%         fprintf(1,'j2=%d, result = %+20.18f  Difference = %+20.18f \n\n',j2,Q2,diff2);
%
%
%     end % for statement

%
% end ZSquare validation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%  Validate ZAnnulus

%     fprintf(1,'\n\n    ANNULUS  VALIDATION \n');
%     for j1=1:37
%
%         j2 = j1 - 1;
%         if j1 == 1
%           j2 = 2;
%         end % if statement
%
%
%         [n1,m1,sc1] = convertj(j1, 'MAHAJAN');
%         [n2,m2,sc2] = convertj(j2, 'MAHAJAN');
%
%         ir = 0.3;  % annulus inner radius
%         fun1 = @(rho,theta) ZAnnulus(j1,n1,m1,rho,theta,ir).*ZAnnulus(j1,n1,m1,rho,theta,ir).*rho;
%         Q1 = quad2d(fun1, ir,1, 0,2*pi);
%         diff1 = (pi - pi*ir^2) - Q1;
%
%         fun2 = @(rho,theta) ZAnnulus(j1,n1,m1,rho,theta,ir).*ZAnnulus(j2,n2,m2,rho,theta,ir).*rho;
%         Q2 = quad2d(fun2, ir,1, 0,2*pi);
%         diff2 = 0.0 - Q2;
%
%         fprintf(1,'j1=%d, n1=%d, m1=%d, result = %+20.18f  Difference = %+20.18f \n',j1,n1,m1,Q1,diff1);
%         fprintf(1,'j2=%d, n2=%d, m2=%d, result = %+20.18f  Difference = %+20.18f \n\n',j2,n2,m2,Q2,diff2);
%
%     end % for statement
%
%
% end annulus validation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end % validateZ
