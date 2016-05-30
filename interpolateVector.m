function [ compx_x_y, compy_x_y ] = interpolateVector( x, y, max_x, max_y, field )
%INTERPOLATEVECTOR This function bilinear interpolate a 2d vector field
%   Return a normalized vector (length = 1), if the original length too
%   small, return vector 0;
%   x, y: coordinate of the desired vector
%   max_x, max_y: resolution of the target space
%   field: the vector field, which will be stretch to match the resolution
%           of the target space.
    %First translate to the field resolution
    function [result] = interpolate(delta, f1, f2)
        % Linear interpolate
        result = f1+(f2-f1)*delta;
    end
    if (x<1 || y<1 ||x>max_x||y>max_y) %Return 0 vector if out of range
        compx_x_y = 0;
        compy_x_y = 0;
        return;
    end
    x = (x-1)/max_x*(size(field, 2)-1)+1;
    y = (y-1)/max_y*(size(field, 3)-1)+1;

    x1 = floor(x);
    x2 = ceil(x);
    y1 = floor(y);
    y2 = ceil(y);
    delta_x = x - x1;
    delta_y = y - y1;
    
    compx_x_y1 = interpolate( delta_x, field(1, x1, y1), field(1, x2, y1));
    compx_x_y2 = interpolate( delta_x, field(1, x1, y2), field(1, x2, y2));
    
    compy_x_y1 = interpolate( delta_x, field(2, x1, y1), field(2, x2, y1));
    compy_x_y2 = interpolate( delta_x, field(2, x1, y2), field(2, x2, y2));
    
    compx_x_y = interpolate( delta_y, compx_x_y1, compx_x_y2);
    compy_x_y = interpolate( delta_y, compy_x_y1, compy_x_y2);
end