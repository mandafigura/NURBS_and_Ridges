function [kappa_1u,kappa_1v,kappa_2u,kappa_2v] = check_critical_condition(x)
    point = the_coeficients_of_a_point_on_surface(x);
    %% Algorithm description
    %- k1 and k2 derivatives to be checked
    kappa_1u = (point.Pu) + ((point.Ru)./sqrt(point.Q));
    kappa_1v = (point.Pv) + ((point.Rv)./sqrt(point.Q));
    kappa_2u = (point.Pu) - ((point.Ru)./sqrt(point.Q));
    kappa_2v = (point.Pv) - ((point.Rv)./sqrt(point.Q));
end