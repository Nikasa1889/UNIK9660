function [ texture_map ] = computeIlluminateTextureMap()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    %% Parameter
    ka = 0.1;
    kd = 0.3;
    ks = 0.6;
    n = 40;
    p = 4.8;
    PRECISION = 1001;
    N = PRECISION;
    %% Compute
    texture_map = zeros(N, N);
    for i = 1:N
        for j = 1:N
            t1 = (i-1)/(N-1);
            t2 = (j-1)/(N-1);
            I_amb = ka;
            LT = 2*t1-1;
            LN = sqrt(1-LT^2);
            I_diff = kd*(LN)^p; %p use for excess brighness problem
            VT = 2*t2-1;
            VR = LT*VT-sqrt(1-LT^2)*sqrt(1-VT^2);
            I_spec = ks*(VR)^n;
            
            texture_map(i, j) = I_amb + I_diff + I_spec;
        end
    end
end