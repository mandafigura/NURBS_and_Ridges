function [DX,DY] = derivative_order_c_NURBS_curve(u,U,P,n,c)
    %% Computes the (c-1)th derivative of a NURBS curve of degree n-1
    % Input:  
    %         Parameter u
    %         Knot Vector U
    %         Defining points P of the curve (homogenous coordinates)
    %         B-Splines degree is n-1
    %         Diff order is c-1
    % Return: 
    %         Derivatives til the order c-1 of the curve in each parameter [DX, DY]
    
    %% Computes all orders derivatives
    for k=1:c
        if k == 1
            %% Computes order 0 derivatives of the curve point C_w(u)
            [DW(1), DW(2), DW(3)] = parameter_u_b_spline_curve(u,U,P,n)
            [H] = H_Perspective(DW);
            % The initial values
            Wu(k) = DW(3);
            DX(k) = H(1)
            DY(k) = H(2)
        else
            %% Computes remaining orders derivatives of the curve point C_w(u)
            % Computes the derivative of all N_i,k
            for i=1:(length(U)-n)
               DN(i,:) =  diff_basis_func_on_u(n,u,U,c,i);
            end
%%
            % Starting values
            DCX(k)=0;
            DCY(k)=0;
            DCZ(k)=0;            
            % Computes derivatives for B-Splines curve
            for j = 1:size(P,2)
                DCX(k) = DCX(k) + DN(j,k)*P(1,j);
                DCY(k) = DCY(k) + DN(j,k)*P(2,j);
                DCZ(k) = DCZ(k) + DN(j,k)*P(3,j);
            end
%%
            % The terms that build C
            Au = [DCX(k), DCY(k)]
            Wu(k) = DCZ(k)
%%
            % Computes the formula sum
            for j = 2:k
                Wu(j)
                Au = Au - (nchoosek(k-1,j-1)*Wu(j)*[DX(k-j+1),DY(k-j+1)]);
            end
%%
            % Derivative formula
            DX(k) = Au(1)/Wu(1);
            DY(k) = Au(2)/Wu(1);
        end
    end
end
