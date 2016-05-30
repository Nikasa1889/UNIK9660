function [ I ] = plotIlluminatedLines( fieldline, texture_map, I, from_idx, max_arclength )
%Plot a redefined fieldline
%This function use HSV color model
    if (~ exist('from_idx', 'var'))
        from_idx = 1;
    end
    if (~ exist('max_arclength', 'var'))
        max_arclength = 100000;
    end
    
    max_x = size(I, 1);
    max_y = size(I, 2);
    [max_t1, max_t2] = size(texture_map);
    height = (max_x+max_y)/2;
    %% Parameters
    L_point = [max_x/2, max_y/2, height];
    V_point = [max_x/2, max_y/2, height];
    Hue = 60/255; Sat = 25/255;
    arclength = 0;
    for i = from_idx:(length(fieldline.x)-1)
        curr_point = [fieldline.x(i), fieldline.y(i), 0];
        next_point = [fieldline.x(i+1), fieldline.y(i+1), 0];
        
        
        T = next_point - curr_point;
        arclength = arclength + sqrt(sum(T.^2));
        T = normalize(T);
        
        
        L = normalize(curr_point-L_point);
        V = normalize(curr_point-V_point);
        
        t1 = 1/2*(L*T'+1)*(max_t1-1)+1;
        t2 = 1/2*(V*T'+1)*(max_t2-1)+1;
        Illu = texture_map(round(t1), round(t2));
        %I = bitmapplot(fieldline.x(i:(i+1)), fieldline.y(i:(i+1)), I, struct('LineWidth',1,'Color',[Hue, Sat, Illu, 1]));
        
        %Matlab secret: using new function for drawing line is much slower
        pline = LinePixels(fieldline.x(i), fieldline.y(i), fieldline.x(i+1), fieldline.y(i+1));
        H_channel = zeros(1, size(pline, 2))+1;
        S_channel = zeros(1, size(pline, 2))+2;
        V_channel = zeros(1, size(pline, 2))+3;
        idx_H = sub2ind(size(I), floor(pline(1,:)), floor(pline(2, :)), H_channel);
        idx_S = sub2ind(size(I), floor(pline(1,:)), floor(pline(2, :)), S_channel);
        idx_V = sub2ind(size(I), floor(pline(1,:)), floor(pline(2, :)), V_channel);
        I(idx_H) = Hue;
        I(idx_S) = Sat;
        I(idx_V) = Illu;
        if (max_arclength <= arclength)
            break;
        end
    end
    
    function v = normalize(v)
        v = v./sqrt(sum(v.^2));
    end
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

