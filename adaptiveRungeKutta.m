function [fieldline] = adaptiveRungeKutta( start_x, start_y, max_x, max_y, step_size, max_length, field, return_eqSpace, arclength_min)
    if (~ exist('return_eqSpace', 'var'))
        return_eqSpace = true; %default will return equal space fieldline
    end
    
    if (~ exist('arclength_min', 'var'))
        arclength_min = 0.1; %Min arclength after each move, prevent too small step near singularity
    end
    
    EPSILON = 0.05;
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
    
    [fieldline_fwd, fwd_length] = adaptiveRungeKuttaOneside(abs(step_size), max_length);
    [fieldline_bwd, bwd_length] = adaptiveRungeKuttaOneside(-abs(step_size), max_length);
    %%Comment the 2 following lines if we don't need equal-space fieldline
    if (return_eqSpace)
        [fieldline_fwd, fwd_length] = interpolateFieldline(step_size, fieldline_fwd, fwd_length);
        [fieldline_bwd, bwd_length] = interpolateFieldline(step_size, fieldline_bwd, bwd_length);
    end
    
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
    
    function [fieldline_eq, n] = interpolateFieldline(step_size, fieldline, n_fieldline)
    %Linear Interpolation at equal arclength for adaptive RungeKutta
        fieldline_eq = struct('x', [], 'y', []);
        max_points = length(fieldline.x);
        fieldline_eq.x = zeros(1, max_points);
        fieldline_eq.y = zeros(1, max_points);
        
        x = fieldline.x(1); y = fieldline.y(1);
        fieldline_eq.x(1) = x;
        fieldline_eq.y(1) = y;
        i = 1;
        n = 1;
        curr_arclength = step_size;
        while i < n_fieldline
           x_h = fieldline.x(i+1);
           y_h = fieldline.y(i+1);
           arclength = sqrt((x_h-x)^2 + (y_h-y)^2);
           if (curr_arclength > arclength)
               curr_arclength = curr_arclength - arclength;
               x = x_h;
               y = y_h;
               i = i + 1;
           else
               h = curr_arclength;
               while ( h <= arclength ) && (n <max_points)
                   alpha = h/arclength;
                   fieldline_eq.x(n+1) = x + (x_h-x)*alpha;
                   fieldline_eq.y(n+1) = y + (y_h-y)*alpha;
                   n = n+1;
                   h = h + step_size;
               end
               curr_arclength = h - arclength;
               x = x_h; y = y_h;
               i = i + 1;
           end
        end
    end
    function [fieldline, length] = adaptiveRungeKuttaOneside(step_size, max_length)
    %RungeKutta marks the out_field with a fieldline in a vector field
    %   start_x, start_y: seed point;
    %   max_x, max_y: resolution of the target space
    %   field: the vector field;
    %   out_field: a field that will be marked by the fieldline
    %   TODO: Detect Singularity, discontinues, and adaptive step size!
    %   Normalize vector length!
        h_max = 50*abs(step_size);
        %h_min = abs(step_size);
        TOL = 0.001; rho = 1;
        h = step_size;
        
        x = start_x; y = start_y;
        fieldline = struct('x', [], 'y', []);
        fieldline.x = zeros(1, ceil(max_length/abs(arclength_min)));
        fieldline.y = zeros(1, ceil(max_length/abs(arclength_min)));
        fieldline.x(1) = x; fieldline.y(1) = y;
        i = 1;
        while (max_length>0 && abs(h) > arclength_min)
            [fx, fy] = interpolateVector(x, y, max_x, max_y, field);
            if (fx == 0 && fy == 0)
                h = h/2;
                continue;
            end
            delta_x1 = h*fx;
            delta_y1 = h*fy;
            if (x+delta_x1 > max_x || y+delta_y1 > max_y)
                h = sign(h)*min(abs((max_x - x)/fx), abs((max_y - y)/fy));
                continue;
            end
            
            [fx_1, fy_1] = interpolateVector(x+delta_x1/2, y+delta_y1/2, max_x, max_y, field);
            if (fx_1 == 0 && fy_1 == 0)
                h = h/2;
                continue;
            end
            delta_x2 = h*fx_1;
            delta_y2 = h*fy_1;
            if (x+delta_x2 > max_x || y+delta_y2 > max_y)
                h = sign(h)*min(abs((max_x - x)/fx_1), abs((max_y - y)/fy_1));
                continue;
            end
            
            [fx_2, fy_2] = interpolateVector(x+delta_x2/2, y+delta_y2/2, max_x, max_y, field);
            if (fx_2 == 0 && fy_2 == 0)
                h = h/2;
                continue;
            end
            delta_x3 = h*fx_2;
            delta_y3 = h*fy_2;
            if (x+delta_x3 > max_x || y+delta_y3 > max_y)
                h = sign(h)*min(abs((max_x - x)/fx_2), abs((max_y - y)/fy_2));
                continue;
            end
            
            [fx_3, fy_3] = interpolateVector(x+delta_x3, y+delta_y3, max_x, max_y, field);
            if (fx_3 == 0 && fy_3 == 0)
                h = h/2;
                continue;
            end

            delta_x4 = h*fx_3;
            delta_y4 = h*fy_3;
            x_h = x + delta_x1/6 + delta_x2/3 + delta_x3/3 + delta_x4/6;
            y_h = y + delta_y1/6 + delta_y2/3 + delta_y3/3 + delta_y4/6;
            if (x_h > max_x || y_h > max_y)
                h = h/2;
                continue;
            end

            
            [fx_h, fy_h] = interpolateVector(x_h, y_h, max_x, max_y, field);
            error_x = (delta_x4-h*fx_h)/6;
            error_y = (delta_y4-h*fy_h)/6;
            
            h_star = min(abs(h*nthroot(abs(rho*TOL/error_x), 5)), abs(h*nthroot(abs(rho*TOL/error_y), 5)));
            if (abs(h) > h_star + EPSILON)
                h = sign(h)*h_star;
                continue;
            else
                arclength = sqrt((x-x_h)^2+(y-y_h)^2);
                if (arclength < arclength_min)
                    break;
                end
                max_length = max_length-arclength;
                x = x_h;
                y = y_h;
                fieldline.x(i+1) = x; fieldline.y(i+1) = y;
                h = sign(h)*min(h_star, h_max);                
                i = i + 1;
                if abs(h) < arclength_min || x > max_x || x < 1 || y > max_y || y < 1 %boundary and singular check
                    break;
                end
            end
        end
        length = i;
    end
end

