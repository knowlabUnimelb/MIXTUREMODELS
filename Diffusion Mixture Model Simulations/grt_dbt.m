% General Recognition Theory and Decision Bound Theory
%  
% Compute drift rate (p(Step to A)) given dimension boundaries, means and variances
% by integrating a bivariate normal distribution (or pair of marginal
% distributions) in the category A region
function [p] = grt_dbt(x, dx, sx, coactive)
% Inputs
%   x - matrix of mean item locations
%       For logical rule studies where other means of estimating the
%       mean location is not available (e.g., MDS)
%       x = [3 3; 3 2; 2 3; 2 2; 3 1; 2 1; 1 3; 1 2; 1 1];
%   dx - decision bound locations (horizontal, vertical)
%       dx = [1.5 1.5]
%   sx - perceptual variability of normal distribution on dimension 1 and
%       dimension 2
%       sx = [.5 .5]
%   coactive - set equal to true to use bivariate normal; set to false to
%       use two marginal distribution (e.g., for serial or parallel models)

if ~coactive
    if size(x, 1) > 1
        dx = repmat(dx, size(x, 1), 1);
        sx = repmat(sx, size(x, 1), 1);
    end
    
    zx = ((dx - x)./sx);
    p  = normcdf(zx, 0, 1);
    
    % Truncate drift rate
    p(p<.00001) = .00001;
    p(p>.99999) = .99999;
    % p2 = ztop(zx);

else
    %   B | D
    %  -------
    %   A | C

   % Check positive semi-definiteness 
   [~,err] = cholcov(diag(sx),0);
   if err~=0
       warning('Covariance matrix not positive definate. Likely cause is that the variances have grown unusually large.');
       save errortemp % Load this .mat to check the problem
   end
    
   % If positive semi-definiteness is ok, set up regions to integrate
    aRegion = [-inf -inf; dx(1) dx(2)];
    bRegion = [-inf dx(2); dx(1)  inf];
    cRegion = [dx(1) -inf;  inf dx(2)];
    dRegion = [dx(1) dx(2);  inf  inf];
    
    for i = 1:size(x, 1)
        p(i,:) = mvncdf(aRegion(1,:), aRegion(2, :), x(i,:), diag(sx)) + mvncdf(bRegion(1,:), bRegion(2, :), x(i,:), diag(sx)) + mvncdf(cRegion(1,:), cRegion(2,:), x(i,:), diag(sx));
%         q(i,:) = mvncdf(dRegion(1,:), dRegion(2, :), x(i,:), diag(sx)); %
%         q = 1 - p
    end 
end
p = 1 - p; 