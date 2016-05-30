function [ I ] = animateStreamLine( max_x, max_y, field )
%This function animate a streamline
    step_size = 0.5;
    max_length = 1000;
    I = zeros([max_x max_y 3]);
    texture_map = computeIlluminateTextureMap();
    fieldlines = struct('x', [], 'y', [], 'fwd_length', 0, 'bwd_length', 0, 'middle_idx', 0);
    q = 0.90;
    figure;
    for i = 1:12
        for j = 1:12
            start_x = i * 50;
            start_y = j * 50;
            fieldline = adaptiveRungeKutta( start_x, start_y, max_x, max_y, step_size, max_length, field, true, 0.05);
            fieldlines((i-1)*12+j) = fieldline;
            %I = plotIlluminatedLines(fieldline, texture_map, I, 1, 100);
        end
    end
    for t = 1:100
        for i = 1:12
            for j = 1:12
                fieldline = fieldlines((i-1)*12+j);
                I = plotIlluminatedLines(fieldline, texture_map, I, (t-1)*10+fieldline.middle_idx, 10);
            end
        end
        imshow(hsv2rgb(I));
        %k = waitforbuttonpress;
        pause(0.05);
        I(:) = q*I(:);
    end
end

