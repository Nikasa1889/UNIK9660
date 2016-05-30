function [ fieldline ] = forwardEuler( start_x, start_y, max_x, max_y, step_size, max_length, field)
%forwardEuler marks the out_field with a fieldline in a vector field
%   start_x, start_y: seed point;
%   field: the vector field;
%   out_field: a field that will be marked by the fieldline
    fieldline = struct('x', [], 'y', [], 'fwd_length', 0, 'bwd_length', 0, 'middle_idx', 0);
    [fx, fy] = interpolateVector(start_x, start_y, max_x, max_y, field);
    if (fx == 0 && fy == 0)
        fieldline.middle_idx = 1;
        fieldline.fwd_length = 1;
        fieldline.bwd_length = 1;
        fieldline.x = start_x;
        fieldline.y = start_y;
        return;
    end
    
    [fieldline_fwd, fwd_length] = forwardEulerOneSide(abs(step_size), max_length);
    [fieldline_bwd, bwd_length] = forwardEulerOneSide(-abs(step_size), max_length);

    fieldline_bwd.x = fliplr(fieldline_bwd.x); %reverse the backward fieldline
    fieldline_bwd.y = fliplr(fieldline_bwd.y);
        
    fieldline.x = [fieldline_bwd.x(1:end-1), fieldline_fwd.x]; %concate forward and backward
    fieldline.y = [fieldline_bwd.y(1:end-1), fieldline_fwd.y];
    %Remove -1
    from_idx = length(fieldline_bwd.x) + 1 - bwd_length;
    to_idx = from_idx + bwd_length+fwd_length-2;
    fieldline.x = fieldline.x(from_idx:to_idx);
    fieldline.y = fieldline.y(from_idx:to_idx);
    
    fieldline.fwd_length = fwd_length;
    fieldline.bwd_length = bwd_length;
    fieldline.middle_idx = bwd_length;
    
    
    function [fieldline, length] = forwardEulerOneSide(step_size, max_length)
        x = start_x; y = start_y;
        fieldline = struct('x', [], 'y', []);
        fieldline.x = zeros(1, ceil(max_length/abs(step_size)));
        fieldline.y = zeros(1, ceil(max_length/abs(step_size)));
        fieldline.x(1) = x; fieldline.y(1) = y;
        i = 1;
        while max_length>0
            [fx, fy] = interpolateVector(x, y, max_x, max_y, field);
            
            if (fx == 0 && fy == 0) %singularity and discontinues check
                break;
            end
            x = x + step_size*fx;
            y = y + step_size*fy;
            
            if (x > max_x || x < 1 || y > max_y || y < 1) %boundary check
                break;
            end
            max_length = max_length - abs(step_size);
            fieldline.x(i+1) = x; fieldline.y(i+1) = y;
            i = i + 1;
        end
        length = i;
    end
end

