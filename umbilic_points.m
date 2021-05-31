function [y] = umbilic_points(x)
    point = the_coeficients_of_a_point_on_surface(x);
    %% Algorithm description
    %--Generates surface umbilic points system of equations
    y(1) = 2.*(point.barB).*(point.barBu) - 4.*(point.Au).*(point.barC) - 4.*(point.A).*(point.barCu);
    y(2) = 2.*(point.barB).*(point.barBv) - 4.*(point.Av).*(point.barC) - 4.*(point.A).*(point.barCv);
end