function [lambda] = baricentric_coordinates_2d(a,b,c,p)
    %% Description:
    %               Finds the baricentric coordinates of a point p \in R^2
    %               with respect to the 2-simplex with coordinates a,b,c.
    %% Output:
    %               lambda is a vector of size 3
    T = [1 1 1; a(1) b(1) c(1); a(2) b(2) c(2)];
    v = [1;p(1);p(2)];
    lambda = T\v;
end