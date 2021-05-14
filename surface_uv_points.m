function [C] = surface_uv_points(n,m,P,u,v,U,V)
    NU = basis_functions_u(u,U,n);
    NV = basis_functions_u(v,V,m);
    C = zeros(1,size(P,3));
    for i = 1:size(P,1)
        for j = 1:size(P,2)
            for h = 1:size(P,3)     %goes through axis
            C(h) = C(h) + NU(i,n).*NV(j,m).*P(i,j,h);
        end
    end        
end