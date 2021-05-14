function [x2,y2,z2,flag2,u2] = check_simplexes2(d11,d12,d21,d22)
    %% Algorithm description
    %- Tries to find a initial simplex for Regula-Falsi
    %% Input:
    %- A convex subset of R^2 (a square with vertices d11,d12,d21,d22)
    %% Output:
    %- The simplex vertices satifying the initial codition for Regula-Falsi
    %- a flag saying if the Regula-Falsi should run (in case its true)
    u2 = [Inf Inf];
    flag2 = false;
    %% Divides the square in two simplexes
    %% Simplex one
    x2 = sis_func2(d11);
    y2 = sis_func2(d12);
    z2 = sis_func2(d22);
    %x = [(d11(1) - 1) (d11(2) -1)];
    %y = [(d12(1) - 1) (d12(2) -1)];
    %z = [(d22(1) - 1) (d22(2) -1)];
    %- Writes p = (0,0) in baricentric coordinates  
    lambda2 = baricentric_coordinates_2d(x2,y2,z2,[0 0]);
    %- Checks conditions to see if (0,0) is in the simplex.
    if(lambda2 >= 0)
        flag2 = true;
        u2 = lambda2(1)*d11 + lambda2(2)*d12 + lambda2(3)*d22;
        x2 = d11;
        y2 = d12;
        z2 = d22;
        return;
    else
        flag2 = false;
    end
    %% Simplex two
    x2 = sis_func2(d11);
    y2 = sis_func2(d21);
    z2 = sis_func2(d22);
    %x = [(d11(1) - 1) (d11(2) -1)];
    %y = [(d21(1) - 1) (d21(2) -1)];
    %z = [(d22(1) - 1) (d22(2) -1)];
    %- Writes p = (0,0) in baricentric coordinates  
    lambda2 = baricentric_coordinates_2d(x2,y2,z2,[0 0]);
    %- Checks conditions to see if (0,0) is in the simplex.
    if(lambda2 >= 0)
        flag2 = true;
        u2 = lambda2(1)*d11 + lambda2(2)*d21 + lambda2(3)*d22;
        x2 = d11;
        y2 = d21;
        z2 = d22;
    else
        flag2 = false;
    end
end