function [D] = derivative_order_cd_bspline_surface(u,v,U,V,P,n,m,c,d)
%% Computes the (c-1)(d-1)th partial derivative of a B-Splines surface of degree (n-1)x(m-1)
    % Input:  
    %         Parameter u
    %         Parameter v
    %         Knot Vector U
    %         Knot Vector V
    %         Defining points P of the curve
    %         B-Splines degrees are (n-1) and (m-1)
    %         Partial-Diff order is (c-1)x(d-1)
    % Return: 
    %         Partial Derivatives of order (c-1)x(d-1) of the surface in each parameter [D]
    
    
    %% Computes de derivatives of each basis functions of degree n and m, for all i,j
    for i=1:(length(U)-n)
        DNU(i,:) =  diff_basis_func_on_u(n,u,U,c,i);
    end
    for j=1:(length(V)-m)
        DNV(j,:) =  diff_basis_func_on_u(m,v,V,d,j);
    end
    %% Computes the derivatives in each axis
    D = zeros(1,size(P,3));
    for i = 1:size(P,1)
        for j = 1:size(P,2)
            for h = 1:size(P,3)     %goes through axis
               D(h) = D(h) + DNU(i,c).*DNV(j,d).*P(i,j,h); 
            end
        end
    end
end