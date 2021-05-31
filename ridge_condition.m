function [ridge_cond_k1,ridge_cond_k2] = ridge_condition(x)
    
    [kappa_1u,kappa_1v,kappa_2u,kappa_2v] = check_critical_condition(x); %principal curvature partial derivatives
    [t1,t2] = principal_directions2d(x);
    
    %blue (maximal) ridge condition
    ridge_cond_k1 = dot([kappa_1u kappa_1v],t1);
    
    %red  (minimal) ridge condition
    ridge_cond_k2 = dot([kappa_2u kappa_2v],t2);
end