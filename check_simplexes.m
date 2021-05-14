function [x,y,z,flag,u] = check_simplexes(d11,d12,d21,d22)
    %% Algorithm description
    %- Tries to find a initial simplex for Regula-Falsi
    %% Input:
    %- A convex subset of R^2 (a square with vertices d11,d12,d21,d22)
    %% Output:
    %- The simplex vertices satifying the initial codition for Regula-Falsi
    %- a flag saying if the Regula-Falsi should run (in case its true)
    u = [Inf Inf];
    flag = false;
    %% Divides the square in two simplexes
    %% Simplex one
    x = sis_func(d11);
    y = sis_func(d12);
    z = sis_func(d22);
    %x = [(d11(1) - 1) (d11(2) -1)];
    %y = [(d12(1) - 1) (d12(2) -1)];
    %z = [(d22(1) - 1) (d22(2) -1)];
    %- Writes p = (0,0) in baricentric coordinates  
    lambda = baricentric_coordinates_2d(x,y,z,[0 0]);
    %- Checks conditions to see if (0,0) is in the simplex.
    if(lambda >= 0)
        flag = true;
        u = lambda(1)*d11 + lambda(2)*d12 + lambda(3)*d22;
        x = d11;
        y = d12;
        z = d22;
        return;
    else
        flag = false;
    end
    %% Simplex two
    x = sis_func(d11);
    y = sis_func(d21);
    z = sis_func(d22);
    %x = [(d11(1) - 1) (d11(2) -1)];
    %y = [(d21(1) - 1) (d21(2) -1)];
    %z = [(d22(1) - 1) (d22(2) -1)];
    %- Writes p = (0,0) in baricentric coordinates  
    lambda = baricentric_coordinates_2d(x,y,z,[0 0]);
    %- Checks conditions to see if (0,0) is in the simplex.
    if(lambda >= 0)
        flag = true;
        u = lambda(1)*d11 + lambda(2)*d21 + lambda(3)*d22;
        x = d11;
        y = d21;
        z = d22;
    else
        flag = false;
    end
end