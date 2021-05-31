function [dy] = critical_points_jacobian(x)
    point = the_coeficients_of_a_point_on_surface(x);
    %% Algorithm description
    %--Generates curvature critical points system of equations Jacobian Matrix
    dy(1,1) = (point.Qu).*((point.Pu).^2) + (point.Q).*2.*(point.Pu).*(point.Puu) - 2.*(point.Ru).*(point.Ruu);
    dy(1,2) = (point.Qv).*((point.Pu).^2) + (point.Q).*2.*(point.Pu).*(point.Puv) - 2.*(point.Ru).*(point.Ruv);
    dy(2,1) = (point.Qu).*((point.Pv).^2) + (point.Q).*2.*(point.Pv).*(point.Pvu) - 2.*(point.Rv).*(point.Rvu);
    dy(2,2) = (point.Qv).*((point.Pv).^2) + (point.Q).*2.*(point.Pv).*(point.Pvv) - 2.*(point.Rv).*(point.Rvv);
end