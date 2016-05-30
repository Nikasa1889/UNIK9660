function [ I ] = plotStreamLinesGrid( max_x, max_y, step_size, field, is_Illuminated, numericalIntergrator)
%This function draws streamlines using grid seed points
%You can choose cho illuminate the streamline or not
%The adaptiveRungeKutta is used when is_Illuminated is true, since it's
%super fast for this task
    I = zeros([max_x max_y 3]);
    max_length = max_x+max_y;
    COLOR = [0 1 1 1];
    if is_Illuminated
        texture_map = computeIlluminateTextureMap();
    end
    space_x = floor(max_x/10);
    space_y = floor(max_y/10);
    for start_x = space_x:space_x:(max_x - space_x)
        for start_y = space_y:space_y:(max_y-space_y)
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
    end
    %Convert HSV to RGB
    figure, imshow(hsv2rgb(I));
end