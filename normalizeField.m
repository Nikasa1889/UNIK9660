function [ field ] = normalizeField( field )
% This function normalize a vector field into unit vectors
    EPSILON = 0.0001;
    arc_length = sqrt((field(1,:,:).^2 + field(2,:,:).^2));
    arc_length(arc_length < EPSILON) = 0;
    field(1,:,:) = field(1,:,:) ./ arc_length;
    field(2,:,:) = field(2,:,:) ./ arc_length;
    field(isnan(field)) = 0;
    field(isinf(field)) = 0;
end

