function [lambda] = baricentric_coordinates_2d(a,b,c,p)
    T = [1 1 1; a(1) b(1) c(1); a(2) b(2) c(2)];
    v = [1;p(1);p(2)];
    lambda = T\v;
end