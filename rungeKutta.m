function [fieldline] = rungeKutta( start_x, start_y, max_x, max_y, step_size, max_length, field)
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
    
    [fieldline_fwd, fwd_length] = rungeKuttaOneside(abs(step_size), max_length);
    [fieldline_bwd, bwd_length] = rungeKuttaOneside(-abs(step_size), max_length);

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
    
    function [fieldline, length] = rungeKuttaOneside(step_size, max_length)
    %RungeKutta marks the out_field with a fieldline in a vector field
    %   start_x, start_y: seed point;
    %   max_x, max_y: resolution of the target space
    %   field: the vector field;
    %   out_field: a field that will be marked by the fieldline
    %   TODO: Detect Singularity, discontinues, and adaptive step size!
        x = start_x; y = start_y;
        fieldline = struct('x', [], 'y', []);
        fieldline.x = zeros(1, ceil(max_length/abs(step_size)));
        fieldline.y = zeros(1, ceil(max_length/abs(step_size)));
        fieldline.x(1) = x; fieldline.y(1) = y;
        i = 1;
        while max_length>0
            [fx, fy] = interpolateVector(x, y, max_x, max_y, field);
            if (fx == 0 && fy == 0)
                break;
            end
            delta_x1 = step_size*fx;
            delta_y1 = step_size*fy;

            [fx_1, fy_1] = interpolateVector(x+delta_x1/2, y+delta_y1/2, max_x, max_y, field);
            if (fx_1 == 0 && fy_1 == 0)
                break;
            end
            delta_x2 = step_size*fx_1;
            delta_y2 = step_size*fy_1;

            [fx_2, fy_2] = interpolateVector(x+delta_x2/2, y+delta_y2/2, max_x, max_y, field);
            if (fx_2 == 0 && fx_2 == 0)
                break;
            end
            delta_x3 = step_size*fx_2;
            delta_y3 = step_size*fy_2;

            [fx_3, fy_3] = interpolateVector(x+delta_x3, y+delta_y3, max_x, max_y, field);
            if (fx_3 == 0 && fx_3 == 0)
                break;
            end
            delta_x4 = step_size*fx_3;
            delta_y4 = step_size*fy_3;
            x_new = x + delta_x1/6 + delta_x2/3 + delta_x3/3 + delta_x4/6;
            y_new = y + delta_y1/6 + delta_y2/3 + delta_y3/3 + delta_y4/6;
            arc_length = sqrt((x_new-x)^2 + (y_new-y)^2);
            max_length = max_length - arc_length;
            x = x_new;
            y = y_new;
            fieldline.x(i+1) = x; fieldline.y(i+1) = y;
            i = i + 1;
            if (x > max_x || x < 1 || y > max_y || y < 1) %boundary check
                break;
            end

        end
        length = i;
    end
end



