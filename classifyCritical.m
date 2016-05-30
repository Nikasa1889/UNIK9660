function [ result ] = classifyCritical( x, y, max_x, max_y, field )
%This function classify critical points in 2D vector field
    result = 0;
    [v1(1), v1(2)] = interpolateVector( x, y, max_x, max_y, field);
    if (v1(1) == 0 && v1(2) == 0) return; end
    [v2(1), v2(2)] = interpolateVector( x+1, y, max_x, max_y, field);
    if (v2(1) == 0 && v2(2) == 0) return; end
    [v3(1), v3(2)] = interpolateVector( x, y+1, max_x, max_y, field);
    a = v2(1) - v1(1);
    d = v3(2) - v1(2);
    b = v3(1) - v1(1);
    c = v2(2) - v1(2);
    
    p = a+d; %Trace
    q = a*d - b*c; %Det
    delta = p^2 - 4*q;
    if (delta>=0)
        if (q > 0) && (p > 0)
            result = 1; %Repelling node (source)
        else
            if (q > 0) && (p < 0)
                result = 2; %Attracting node (sink)
            else
                result = 3; %Saddle point %There is one example
            end
        end
    else
        if (p == 0)
            result = 4; %Center
        else
            if (p > 0)
                result = 5; %Repelling spiral (source)
            else
                result = 6; %Attracting spiral (sink)
            end
        end
    end
end

