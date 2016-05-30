function [ I ] = plotVectorField( field, I, space_x, space_y, vector_length )
% Plot vectors on a specific grid
% The color of vector is always blue
    max_x = size(I, 1);
    max_y = size(I, 2);
    [X, Y] = meshgrid((vector_length+1):space_x:(max_x-vector_length), (vector_length+1):space_y:(max_y-vector_length));
    i = 1;
    while (i <= length(X(:)))
        x1 = X(i);
        y1 = Y(i);
        [fx, fy] = interpolateVector(x1, y1, max_x, max_y, field);
        x2 = x1+fx*vector_length;
        y2 = y1+fy*vector_length;
        pline = LinePixels(x1, y1, x2, y2);
        
        blue_channel = repmat(3, 1, length(pline(1,:)));
        idx1 = sub2ind(size(I), ceil(pline(1,:)), floor(pline(2, :)), blue_channel);
        idx2 = sub2ind(size(I), floor(pline(1,:)), ceil(pline(2, :)), blue_channel);
        idx3 = sub2ind(size(I), floor(pline(1,:)), floor(pline(2, :)), blue_channel);
        idx4 = sub2ind(size(I), ceil(pline(1,:)), ceil(pline(2, :)), blue_channel);

        I([idx1 idx2 idx3 idx4]) = 1; 
        i = i+1;
    end
    I = rgb2hsv(I);
    figure, imshow(hsv2rgb(I));
    function pline = LinePixels(x1,y1, x2, y2)
    % Calculate the pixels needed to construct a line of 1 pixel thickness
    % between two coordinates.
        xp=[x1 x2];  yp=[y1 y2]; 
        dx=abs(xp(2)-xp(1)); dy=abs(yp(2)-yp(1));
        if(dx==dy)
         if(xp(2)>xp(1)), xline=xp(1):xp(2); else xline=xp(1):-1:xp(2); end
         if(yp(2)>yp(1)), yline=yp(1):yp(2); else yline=yp(1):-1:yp(2); end
        elseif(dx>dy)
         if(xp(2)>xp(1)), xline=xp(1):xp(2); else xline=xp(1):-1:xp(2); end
         yline=linspace(yp(1),yp(2),length(xline));
        else
         if(yp(2)>yp(1)), yline=yp(1):yp(2); else yline=yp(1):-1:yp(2); end
         xline=linspace(xp(1),xp(2),length(yline));   
        end
        pline(1,:)=xline;
        pline(2,:)=yline;
    end
end