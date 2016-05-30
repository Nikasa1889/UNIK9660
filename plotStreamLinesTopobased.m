function [ I ] = plotStreamLinesTopobased(  TopoDegree, max_x, max_y, step_size, field, is_Illuminated, numericalIntergrator )
%Topology-based, for simplicity, I only implement for saddle and spinral
%point, with fixed direction.
    [xs, ys] = find(TopoDegree == 1); %find saddle point
    n_seeds = 15;
    dist = 10;
    n = 0;
    for i = 1:length(xs)
        x = xs(i);
        y = ys(i);
        for j = 1:n_seeds
            Start_x(n+1) = x+dist*j;
            Start_y(n+1) = y;
            Start_x(n+2) = x;
            Start_y(n+2) = y+dist*j;
            Start_x(n+3) = x-dist*j;
            Start_y(n+3) = y;
            Start_x(n+4) = x;
            Start_y(n+4) = y-dist*j;
            n = n + 4;
        end
    end
    [xs, ys] = find(TopoDegree == 2); %find spiral point
    for i = 1:length(xs)
        x = xs(i);
        y = ys(i);
        for j = 1:n_seeds
            Start_x(n+1) = x+dist*j;
            Start_y(n+1) = y+dist*j;
            n = n + 1;
        end
    end
    
    I = zeros([max_x max_y 3]);
    max_length = max_x+max_y;
    COLOR = [90/255 1 1 1];

    if is_Illuminated
        texture_map = computeIlluminateTextureMap();
    end
    for i = 1:n
        start_x = Start_x(i);
        start_y = Start_y(i);
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
end

