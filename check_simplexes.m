function [x,y,z,flag] = check_simplexes(the_function,d11,d12,d21,d22)
    %% Algorithm description
    %- Tries to find a initial simplex for Regula-Falsi
    %% Input:
    %- The function sis_func to check if a root is inside the simplex
    %- A convex subset of R^2 (a square with vertices d11,d12,d21,d22):
    %% Output:
    %- The simplex vertices satifying the initial codition for Regula-Falsi
    %- a flag saying if the Regula-Falsi method should run (true) or not (false)
    
    %% Initial value
    flag = false;
    %% Divides the square in two simplexes
    %% The simplexes to be checked illustration:
    %
    %                d11-----d12        \      d11-----d12
    %                 |       |      ----\      |  \  1 |     (the two
    %                 |       |      ----/      | 2  \  |      simplexes)
    %                d21-----d22        /      d21-----d22
    %    
    %% Simplex one
    x = the_function(d11);
    y = the_function(d12);
    z = the_function(d22);
    %- Writes p = (0,0) in baricentric coordinates 
    lambda = baricentric_coordinates_2d(x,y,z,[0 0]);
    %- Checks conditions to see if (0,0) is in the simplex.
    if(lambda >= 0)
        flag = true;
        x = d11;
        y = d12;
        z = d22;
        return
    end
    %% Simplex two
    x = the_function(d11);
    y = the_function(d21);
    z = the_function(d22);
    %- Writes p = (0,0) in baricentric coordinates  
    lambda = baricentric_coordinates_2d(x,y,z,[0 0]);
    %- Checks conditions to see if (0,0) is in the simplex.
    if(lambda >= 0)
        flag = true;
        x = d11;
        y = d21;
        z = d22;
        return
    end
end