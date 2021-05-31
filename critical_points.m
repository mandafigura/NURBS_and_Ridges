function [y] = critical_points(x)
    point = the_coeficients_of_a_point_on_surface(x);
    %% Algorithm description
    %--Generates curvature critical points system of equations
    y(1) = (point.Q).*((point.Pu).^2) - ((point.Ru).^2);
    y(2) = (point.Q).*((point.Pv).^2) - ((point.Rv).^2);
end