function [CX, CY, CZ] = parameter_u_b_spline_curve(u,U,P,n)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%Computes b-splines curve point for parameter u
    %%Input:  
    %%        Parameter u
    %%        Knot Vector U
    %%        Defining points P of the curve (homogenous coordinates)
    %%        B-Splines degree is n-1
    %%Return: 
    %%        Point of the curve evaluated in u
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Starting values
    NU = basis_functions_u(u,U,n);   %basis functions of degree n-1
    CX = 0;
    CY = 0;
    CZ = 0;
    %The curve definition
    for i = 1:size(P,2)              %computes for number of control points
        CX = CX + NU(i,n).*P(1,i);
        CY = CY + NU(i,n).*P(2,i);
        CZ = CZ + NU(i,n).*P(3,i);
    end        
end