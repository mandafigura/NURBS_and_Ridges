function [dy] = umbilic_points_jacobian(x)
    point = the_coeficients_of_a_point_on_surface(x);
    %% Algorithm description
    %--Generates surface umbilic points system of equations Jacobian Matrix
    dy(1,1) = 2.*(point.barBu).*(point.barBu) + 2.*(point.barB).*(point.barBuu) - 4.*(point.Auu).*(point.barC) - 4.*(point.Au).*(point.barCu) - 4.*(point.Au).*(point.barCu) - 4.*(point.A).*(point.barCuu);
    dy(1,2) = 2.*(point.barBv).*(point.barBu) + 2.*(point.barB).*(point.barBuv) - 4.*(point.Auv).*(point.barC) - 4.*(point.Au).*(point.barCv) - 4.*(point.Av).*(point.barCu) - 4.*(point.A).*(point.barCuv);
    dy(2,1) = 2.*(point.barBv).*(point.barBu) + 2.*(point.barB).*(point.barBuv) - 4.*(point.Auv).*(point.barC) - 4.*(point.Au).*(point.barCv) - 4.*(point.Av).*(point.barCu) - 4.*(point.A).*(point.barCuv);
    dy(2,2) = 2.*(point.barBv).*(point.barBv) + 2.*(point.barB).*(point.barBvv) - 4.*(point.Avv).*(point.barC) - 4.*(point.Av).*(point.barCv) - 4.*(point.Av).*(point.barCv) - 4.*(point.A).*(point.barCvv);
end