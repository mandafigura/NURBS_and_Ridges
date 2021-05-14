function [DS] = derivative_order_cd_NURBS_surface_test(u,v,U,V,P,n,m,c,d)
    %% Computes the (c-1)(d-1)th partial derivative of a NURBS surface of degree (n-1)x(m-1)
    % Input:  
    %         Parameter u
    %         Parameter v
    %         Knot Vector U
    %         Knot Vector V
    %         Defining points P of the curve (homogenous coordinates)
    %         B-Splines degrees are (n-1) and (m-1)
    %         Partial-Diff order is (c-1)x(d-1)
    % Return: 
    %         Partial Derivatives of order (c-1)x(d-1) of the surface in each parameter [D]
    A  = zeros(c+1,d+1);
    W  = zeros(c+1,d+1);
    DS = zeros(c+1,d+1);
    for k=1:c
        for l=1:d
            aux = derivative_order_cd_bspline_surface(u,v,U,V,P,n,m,k,l);
            A(k,l,1:3) = aux(1:3);
            W(k,l,1) = aux(4);
        end
    end
    
    %% order (0,0) derivative (original func)
    aux = surface_uv_points(n,m,P,u,v,U,V);
    H = H_Perspective(aux);
    DS(1,1,1:3) = H;
    
    D(1)  = DS(1,1,1);
    D(2)  = DS(1,1,2);
    D(3)  = DS(1,1,3);
    if(c==1 & d==1)
        return;
    end
   
    %% order (1,0) derivative
    DS(2,1,1:3) = (A(2,1,1:3) - (W(2,1,1).*DS(1,1,1:3)))./W(1,1,1);
    D(1) = DS(2,1,1);
    D(2) = DS(2,1,2);
    D(3) = DS(2,1,3);
    if(c==2 & d==1)
        return;
    end
    
    %% order (0,1) derivative
    DS(1,2,1:3) = (A(1,2,1:3) - (W(1,2,1).*DS(1,1,1:3)))./W(1,1,1);
    D(1) = DS(1,2,1);
    D(2) = DS(1,2,2);
    D(3) = DS(1,2,3);
    if(c==1 & d==2)
        return;
    end
    
    
       %% derivative of order (k,l)
        for k=1:c
            for l=1:d
                left_sum = 0;
                for i=2:k
                    left_sum = left_sum + nchoosek(k-1,i-1).*W(i,1,1).*DS(k-i+1,l,1:3);
                end
                right_sum = 0;
                for j=2:l
                    right_sum = right_sum + nchoosek(l-1,j-1).*W(1,j,1).*DS(k,l-j+1,1:3);
                end
                middle_sum = 0;
                for i=2:k
                    middle_sum = nchoosek(k-1,i-1).*middle_sum;
                    for j=2:l
                        middle_sum = middle_sum + nchoosek(l-1,j-1).*W(i,j,1).*DS(k-i+1,l-j+1,1:3);
                    end
                end
                DS(k,l,1:3) = (A(k,l,1:3) - left_sum - middle_sum - right_sum)./W(1,1,1);
            end
        end
    D(1) = DS(c,d,1);
    D(2) = DS(c,d,2);
    D(3) = DS(c,d,3);
end