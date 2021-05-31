function [discriminant] = check_umbilic_condition(x)
    point = the_coeficients_of_a_point_on_surface(x);
    %% Algorithm description
    %- k1 = k2 (umbilic) condition to be checked
    discriminant = (point.barB).^2 - 4.*(point.A).*(point.barC);
end