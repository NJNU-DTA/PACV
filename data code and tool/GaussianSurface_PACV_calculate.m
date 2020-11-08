%% Calculate Gaussian surface and its partial derivatives
A = 3;
B = 10;
C = 1/3;
m = 500;
n = 500;
xi = -500:5:500;
yi = -500:5:500;
[x, y] =meshgrid(xi, yi); 

% This can generate a DEM data with 5m resolution.
z = A.*(1-((x./m).^2)).*(exp(-((x./m).^2)-(((y./n)+1).^2))) ...
    - B.*((0.2).*(x./m)-((x./m).^3)-((y./n).^5)).*(exp(-((x./m).^2)-((y./n).^2))) ...
    - C.*(exp(-(((x./m)+1).^2)-((y./n).^2)));


% The partial derivative of Gaussian Surface with respect to X.
fx = exp(- x.^2./250000 - y.^2./250000).*((3.*x.^2)./12500000 - 1./250) ...
    + (exp(- (x./500 + 1).^2 - y.^2./250000).*(x./125000 + 1./250))./3 ...
    - (3.*x.*exp(- (y./500 + 1).^2 - x.^2./250000))./125000 ...
    - (x.*exp(- x.^2./250000 - y.^2./250000).*(x.^3./12500000 - x./250 ...
    + y.^5./3125000000000))./125000 + (x.*exp(- (y./500 + 1).^2 ...
    - x.^2./250000).*((3.*x.^2)./250000 - 3))./125000;


% The partial derivative of Gaussian Surface with respect to Y.
fy =(y.^4.*exp(- x.^2./250000 - y.^2./250000))./625000000000 ...
    + (y.*exp(- (x./500 + 1).^2 - y.^2./250000))./375000 ...
    + exp(- (y./500 + 1).^2 - x.^2./250000).*(y./125000 ...
    + 1./250).*((3.*x.^2)./250000 - 3) - (y.*exp(- x.^2./250000 ...
    - y.^2./250000).*(x.^3./12500000 - x./250 + y.^5./3125000000000))./125000;
 

%% Calculate aspect
AspectMatrix = 57.29578 .* atan2(fy, -fx);
IndexLess0 = find(AspectMatrix < 0);
IndexGreat90 = find(AspectMatrix > 90);
IndexElse = find((AspectMatrix>=0) & (AspectMatrix <=90));
 for i = IndexLess0
     AspectMatrix(i) = 90 - AspectMatrix(i);
 end
 for i = IndexGreat90
     AspectMatrix(i) = 450 - AspectMatrix(i);
 end
 for i = IndexElse
     AspectMatrix(i) = 90 - AspectMatrix(i);
 end
 
  
%% algorothm process
Theta = AspectMatrix;
Alpha = zeros(size(Theta));

% see the formula (1) in paper.
Index_1 = find(Theta >= 0 & Theta <= 90);
 for i = Index_1
     Alpha(i) = 90 - Theta(i);
 end
Index_2 = find(Theta > 90 & Theta <= 360);
for i = Index_2
     Alpha(i) = 450 - Theta(i);
end

cellsize = 5;
 

% see the formula (2) in paper.
Vector_x = cellsize.*cos(Alpha.*(pi./180));
Vector_y = cellsize.*sin(Alpha.*(pi./180));
Vector_asp = cat(3, Vector_x, Vector_y);

[m, n] = size(AspectMatrix);
PACVmatrix = zeros(m, n);

%Extend the boundary of the matrix by replicating.
Vector_asp = padarray(Vector_asp,[1 1], 'replicate');

% slide window
 for i = 2 : m + 1
     for j = 2 : n + 1
         tempValue = FiniteDifference(i, j , Vector_asp, cellsize);
         PACVmatrix(i - 1, j - 1) = tempValue;
     end
 end  

%display the PACV result
imshow(PACVmatrix,[])
 
 %% The FiniteDifference function
 function [out] = FiniteDifference(cr, cc, vm, cellsize)
%FINITEDIFFERENCE the third-order finite difference method (also called Horn¡¯s method)
%   cr is the row number of the central cell.
%   cc is the column number of the central cell.
%   vm is the vector matrix that consist of (X, Y)
%   cellsize is the resolution of the raster.

    p = zeros(8, 2);
    
        
    % we contrust a local 3*3 window and give a index like below; 
    %    6|7|8
    %    5|*|1
    %    4|3|2   
    p(:, 1) = [vm(cr, cc + 1, 1) vm(cr + 1, cc + 1, 1) vm(cr + 1, cc, 1)...
        vm(cr + 1, cc - 1, 1) vm(cr, cc - 1, 1) vm(cr - 1, cc - 1, 1)...
        vm(cr - 1, cc, 1) vm(cr - 1, cc + 1, 1)];
    p(:, 2) = [vm(cr, cc + 1, 2) vm(cr + 1, cc + 1, 2) vm(cr + 1, cc, 2)...
        vm(cr + 1, cc - 1, 2) vm(cr, cc - 1, 2) vm(cr - 1, cc - 1, 2)...
        vm(cr - 1, cc, 2) vm(cr - 1, cc + 1, 2)];

    tPara = 8.*cellsize;
    
    % we contrust a local 3*3 window and give a index like below; 
    %    6|7|8
    %    5|*|1
    %    4|3|2
    
    %see Formulas (3) and (4).
    ACV_we = ((p(8,:)-p(6,:)) + 2*(p(1,:)-p(5,:)) +(p(2,:)-p(4,:)))./tPara;
    ACV_ns = ((p(8,:)-p(2,:)) + 2*(p(7,:)-p(3,:)) +(p(6,:)-p(4,:)))./tPara;

    %see Formula (5) and (6).
    out = ACV_we(1) + ACV_ns(2);
 end

