function [ I ] = plotStreamLinesProp( p, max_x, max_y, step_size, field, is_Illuminated, numericalIntergrator )
%This function is similar to plotStreamLinesGrid, but generate seed points
%based on property p (degree of interest), using Monte Carlo process.
    n_seeds = 100;
    I = zeros([max_x max_y 3]);
    max_length = max_x+max_y;
    p_acc = zeros(1, max_x*max_y);
    COLOR = [60/255 1 1 1];
    if is_Illuminated
        texture_map = computeIlluminateTextureMap();
    end

    p_acc(1) = p(1);
    for i = 2:max_x*max_y
        p_acc(i) = p_acc(i-1) + p(i);
    end
    p_max = p_acc(max_x*max_y);
    for i = 1:n_seeds
        k = unifrnd(0, p_max);
        [v, idx] = closest_value(p_acc, k); 
        [start_x, start_y] = ind2sub([max_x max_y], idx);
        if is_Illuminated
                fieldline = adaptiveRungeKutta( start_x, start_y, max_x, max_y, step_size, max_length, field, false, 0.1);
                I = plotIlluminatedLines(fieldline, texture_map, I);
            else
                %If using adaptiveRungeKutta, we return non-equal space fieldline
                if strcmp(func2str(numericalIntergrator), 'adaptiveRungeKutta') 
                    fieldline = adaptiveRungeKutta( start_x, start_y, max_x, max_y, step_size, max_length, field, false, 0.1);
                else
                    fieldline = numericalIntergrator( start_x, start_y, max_x, max_y, step_size, max_length, field);
                end
                I = bitmapplot(fieldline.x, fieldline.y, I, struct('LineWidth',1,'Color',COLOR));
        end
    end
    figure, imshow(hsv2rgb(I));

    function [v, inf] = closest_value(arr, val)
    % Returns value and index of arr that is closest to val. If several entries
    % are equally close, return the first. Works fine up to machine error (e.g.
    % [v, i] = closest_value([4.8, 5], 4.9) will return [5, 2], since in float
    % representation 4.9 is strictly closer to 5 than 4.8).
    % ===============
    % Parameter list:
    % ===============
    % arr : increasingly ordered array
    % val : scalar in R
        len = length(arr);
        inf = 1;
        sup = len;
        % Binary search for index
        while sup - inf > 1
            med = floor((sup + inf)/2);
            % Replace >= here with > to obtain the last index instead of the first.
            if arr(med) >= val 
                sup = med;
            else
                inf = med;
            end
        end
        % Replace < here with <= to obtain the last index instead of the first.
        if sup - inf == 1 && abs(arr(sup) - val) < abs(arr(inf) - val)
            inf = sup;
        end  

        v = arr(inf);
    end
end

