function [C] = surface_points(n,m,P,u,v,U,V)
%% Computes the (c-1)(d-1)th partial derivative of a B-Splines surface of degree (n-1)x(m-1)
    % Input:  
    %         Parameter u
    %         Parameter v
    %         Knot Vector U
    %         Knot Vector V
    %         Defining points P of the curve (homogenous coordinates)
    %         B-Splines degrees are (n-1) and (m-1)
    % Return: 
    %         The point of the surface on the parameter (u,v)
    
    NU = basis_func(u,U,n);
    NV = basis_func(v,V,m);
    C = zeros(1,size(P,3));
    for i = 1:size(P,1)
        for j = 1:size(P,2)
            for h = 1:size(P,3)    %goes through axis
                C(h) = C(h) + NU(i,n).*NV(j,m).*P(i,j,h);
            end
        end
    end        
end