function [k1,k2] = principal_curvatures(x)
    point = the_coeficients_of_a_point_on_surface(x);
    
    k1 = (-(point.B) + sqrt((point.B).^2 - 4.*(point.A).*(point.C)))./(2.*(point.A));
    k2 = (-(point.B) - sqrt((point.B).^2 - 4.*(point.A).*(point.C)))./(2.*(point.A));
end