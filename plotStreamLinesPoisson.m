function [ I ] = plotStreamLinesPoisson( max_x, max_y, step_size, field, is_Illuminated, numericalIntergrator )
% This function use poisson disc distribution to generate evenly spaced
% seed ponints
    I = zeros([max_x max_y 3]);
    max_length = max_x+max_y;
    min_dist = (max_x + max_y)/20;
    COLOR = [30/255 1 1 1];
    if is_Illuminated
        texture_map = computeIlluminateTextureMap();
    end
    seed_points =  generate_poisson_2d([max_x, max_y], min_dist, 10);
    for i = 1:size(seed_points, 1)
        start_x = floor(seed_points(i, 1));
        start_y = floor(seed_points(i, 2));
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
    %Convert HSV to RGB
    figure, imshow(hsv2rgb(I));

end

