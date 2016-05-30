function [ output ] = lineInterConv( input_texture, field, numericalIntegrator)
%This is a mofified version of Fast LIC algorithm
    tic;
    %% Constants
    [max_x, max_y] = size(input_texture);
    min_Hits = 1;
    min_Points = 5;
    step_size = 0.5;
    L = round(max(max_x/20, max_y/20));
    %L = 15;
    M = 100;
    covered_Threshold = 0.9 * max_x * max_y;
    
    %%
    function [I] = convolute(fieldline, input_texture, from_idx, to_idx)
        k_points = to_idx - from_idx + 1;
        xs = floor(fieldline.x(from_idx:to_idx));
        ys = floor(fieldline.y(from_idx:to_idx));
        ids = sub2ind(size(input_texture), xs, ys);
        I = sum(input_texture(ids));
        I = I / k_points;
    end
    output = struct('I', [], 'numHits', []);
    output.I = zeros(max_x-1, max_y-1);
    output.numHits = zeros(max_x-1, max_y-1);
    %Split image into 10x10 blocks
    n_blockx = 25; n_blocky = 25;
    n_block = n_blockx*n_blocky;
    block_sizex = max_x/n_blockx;
    block_sizey = max_y/n_blocky;
    i = 1;
    n_covered = 0;
    %Error i = 9779
    while i < max_x*max_y
        block_id = mod(i, n_block);
        block_coord = floor(i / (n_block));
        
        block_x = mod(block_id, n_blockx);
        block_y = floor(block_id/(n_blockx));
        
        block_coordx = mod(block_coord, block_sizex);
        block_coordy = floor(block_coord/(block_sizex));
        
        x = (block_x)*block_sizex + block_coordx + 1;
        y = (block_y)*block_sizey + block_coordy + 1;
        if (x~=max_x && y~=max_y)
            if (output.numHits(x, y) >= min_Hits) 
                i = i + 1;
                continue; 
            end
            %Count number of hits, when it's bigger than 90%, use exact convolute
            %Small bug in adaptiveRungeKutta!
            %[fieldline] = forwardEuler(x+0.5, y+0.5, max_x, max_y, step_size, (L+M+1)*step_size, field);
            [fieldline] = numericalIntegrator(x+0.5, y+0.5, max_x, max_y, step_size, (L+M+1)*step_size, field);
            if (length(fieldline.x) <= min_Points)
                i = i + 1;
                continue;
            end
            fwd_length = fieldline.fwd_length;
            bwd_length = fieldline.bwd_length;
            middle_idx = fieldline.middle_idx;
            if fwd_length > L && bwd_length > L && n_covered < covered_Threshold
                max_m_bwd = min(M, bwd_length-L-1);
                max_m_fwd = min(M, fwd_length-L-1);
                middle_idx = middle_idx - max_m_bwd;
                from_idx = middle_idx - L;
                to_idx = middle_idx + L;
            else
                max_m_bwd = 0;
                max_m_fwd = 0;
                from_idx = max(middle_idx - L, middle_idx - bwd_length+1);
                to_idx = min(middle_idx + L, middle_idx + fwd_length-1);
            end
            I0 = convolute(fieldline, input_texture, from_idx, to_idx);
            output.I(x, y) = I0;
            n_covered = n_covered + 1;
            output.numHits(x,y) = output.numHits(x, y) + 1;
            
            k = (to_idx - from_idx);

            for m = (-max_m_bwd+1):max_m_fwd %This loop only happen when fieldline long enough
                x_old = fieldline.x(from_idx);
                y_old = fieldline.y(from_idx);
                from_idx = from_idx + 1;
                to_idx = to_idx + 1;
                middle_idx = middle_idx + 1;
                x = floor(fieldline.x(middle_idx));
                y = floor(fieldline.y(middle_idx));
                x_new = fieldline.x(to_idx);
                y_new = fieldline.y(to_idx);
                
                I0 = I0 + (input_texture(floor(x_new), floor(y_new)) - input_texture(floor(x_old), floor(y_old)))/k;
                if output.numHits(x, y) == 0
                    n_covered = n_covered + 1;
                    if (n_covered >= covered_Threshold)
                        M = 0;
                    end;
                end
                output.I(x, y) = output.I(x, y) + I0;
                output.numHits(x, y) = output.numHits(x, y) + 1;
            end
        end
        i = i + 1;
        if (mod(i, 1000)==0)
            disp(i);
        end
    end
    output.I = output.I ./ output.numHits;
    toc;
end

