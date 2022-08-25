% ***********************************************************************************
% Normalize2dpts - normalizes 2D homogeneous points
%
% Function translates and normalises a set of 2D homogeneous points so that their centroid is at the origin 
% and their mean distance from the origin is sqrt(2).  
% This process typically improves the conditioning of any equations used to solve homographies, fundamental matrices etc.
%
% Usage:   [newpts, T] = normalize2dpts(pts)
% Argument:
%   pts -  3xN array of 2D homogeneous coordinates
% Returns:
%   newpts -  3xN array of transformed 2D homogeneous coordinates.  The
%             scaling parameter is normalised to 1 unless the point is at
%             infinity. 
%   T      -  The 3x3 transformation matrix, newpts = T*pts
%      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [newpts, T] = normalize2dpts(pts,c)

    if size(pts,1) ~= 3
        error('pts must be 3xN');
    end
    
    % Find the indices of the points that are not at infinity
    finiteind = find(abs(pts(3,:)) > eps);
    
    if length(finiteind) ~= size(pts,2)
        warning('Some points are at infinity');
    end
    
    % For the finite points ensure homogeneous coords have scale of 1
    pts(1,finiteind) = pts(1,finiteind)./pts(3,finiteind);
    pts(2,finiteind) = pts(2,finiteind)./pts(3,finiteind);
    pts(3,finiteind) = 1;
    if nargin == 1
    c = mean(pts(1:2,finiteind)')';            % Centroid of finite points;
    end 
    
    newp(1,finiteind) = pts(1,finiteind)-c(1); % Shift origin to centroid.
    newp(2,finiteind) = pts(2,finiteind)-c(2);
    
    meandist = mean(sqrt(newp(1,finiteind).^2 + newp(2,finiteind).^2));
    
    scale = sqrt(2)/meandist;
    
    T = [scale   0   -scale*c(1)
         0     scale -scale*c(2)
         0       0      1      ];
    
    newpts = T*pts;
    
    
    