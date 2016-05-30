function [ I ] = calProperty( I, field, property_name )
%This function calculate a particular property of a 2d vector field
%property_name could be:    "magnitude_Velocity"
%                           "Q_criteria"
%                           "enstrophy"
%                           "topoDegree"
    max_x = size(I, 1);
    max_y = size(I, 2);
    x = 1:max_x+0.5;
    y = 1:max_y+0.5;
    [X, Y]= meshgrid(x, y);
    switch property_name
        case 'magnitude_Velocity'
            output = arrayfun(@magnitude_Velocity, X, Y);
        case 'Q_criteria'
            output = arrayfun(@Q_criteria, X, Y);
        case 'enstrophy'
            output = arrayfun(@enstrophy, X, Y);
        case 'topoDegree'
            output = arrayfun(@topoDegree, X, Y);
        otherwise
            disp('property_name must be "magnitude_Velocity", "Q_criteria", "enstrophy" or "topoDegree"');
            return;
    end
    I = reshape(output, max_x, max_y)';
    return;
    
    function [ result ] = magnitude_Velocity( x, y)
     %This function calculate magnitude of velocity at a given point
        [v_x, v_y] = interpolateVector( x, y, max_x, max_y, field);
        result = sqrt(v_x^2 + v_y^2);
    end
    function [ result ] = Q_criteria(x, y)
     %This function calculate the Q criteria for 2D vector field
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
        Q = -a^2 - d^2 - 2*b*c;
        if Q > 0
            result = Q;
        end
    end

    function [ result ] = enstrophy( x, y )
    %This function calculate an enstrophy of 2D vector field
        result = 0;
        [v1(1), v1(2)] = interpolateVector( x, y, max_x, max_y, field);
        if (v1(1) == 0 && v1(2) == 0) return; end
        [v2(1), v2(2)] = interpolateVector( x+1, y, max_x, max_y, field);
        if (v2(1) == 0 && v2(2) == 0) return; end
        [v3(1), v3(2)] = interpolateVector( x, y+1, max_x, max_y, field);
        dvx_dy = v3(1)-v1(1);
        dvy_dx = v2(2)-v1(2);
        result = (dvx_dy-dvy_dx)^2;
    end
    function [ degree ] = topoDegree( x, y )
    %This function calculate the topological degree of a voxel at x, y
    %The width and length of voxel are fixed at one.
    %x, y must be represented as a middle point of the grid
        EPSILON = 0.05;
        degree = 0;
        %Calculate topological degree
        [v1(1), v1(2)] = interpolateVector( floor(x)-1, floor(y)-1, max_x, max_y, field);
        if (v1(1) == 0 && v1(2) == 0) return; end
        [v2(1), v2(2)] = interpolateVector( floor(x)-1, ceil(y)+1, max_x, max_y, field);
        if (v2(1) == 0 && v2(2) == 0) return; end
        [v3(1), v3(2)] = interpolateVector( ceil(x)+1, ceil(y)+1, max_x, max_y, field);
        if (v3(1) == 0 && v3(2) == 0) return; end
        [v4(1), v4(2)] = interpolateVector( ceil(x)+1, floor(y)-1, max_x, max_y, field);
        if (v4(1) == 0 && v4(2) == 0) return; end  
        theta1 = acos(dot(v1, v2));
        theta2 = acos(dot(v2, v3));
        theta3 = acos(dot(v3, v4));
        theta4 = acos(dot(v4, v1));

        n = (theta1+theta2+theta3+theta4)/(2*pi);
        if abs(round(n) - n) < EPSILON
            degree = round(n);
        else
            degree = 0;
        end
    end
end

