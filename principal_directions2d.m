function [t1,t2] = principal_directions2d(x)
    point   = the_coeficients_of_a_point_on_surface(x);
    [k1,k2] = principal_curvatures(x);
    
    t11 = k1.*(point.F) - (point.f);
    t21 = (point.e) - k1.*(point.E);
    t12 = k2.*(point.F) - (point.f);
    t22 = (point.e) - k2.*(point.E);
    t1op1 = [t11 t21];
    t2op1 = [t12 t22];
    %optional
    t11 = k1.*(point.G) - (point.g);
    t21 = (point.f) - k1.*(point.F);
    t12 = k2.*(point.G) - (point.g);
    t22 = (point.f) - k2.*(point.F);
    t1op2 = [t11 t21];
    t2op2 = [t12 t22];
    
    %which one is better?
    if cond(t1op1) <= cond(t1op2)
        t1 = t1op1;
    else
        t1 = t1op2;
    end

    if cond(t2op1) <= cond(t2op2)
        t2 = t2op1;
    else
        t2 = t2op2;
    end
end