function [x,y,z] = check_simplexes(d11,d12,d21,d22)
    %% Algorithm description
    %- Generates a set of initial points for the Regula-Falsi algorithm;
    %% Input:
    %- A convex subset of R^2 (a square with vertices d11,d12,d21,d22)
    %% Output:
    %- The simplex vertices satifying the initial codition for Regula-Falsi
    
    %% Tries to find a initial simplex for Regula-Falsi
    %- divides the square in two simplexes:
    %  Simplex one
    f(:,1) = sis_func(d11)
    f(:,2) = sis_func(d12)
    f(:,3) = sis_func(d22)
    %- Writes p = (0,0) in baricentric coordinates  
    T = [(f(:,1)-f(:,3)) (f(:,2)-f(:,3))]
    lambda = inv(T)*(-(f(:,3)))
    %- Checks conditions to see if (0,0) is in the simplex.
    if(lambda >= 0 & sum(lambda)<=1)
        regula_falsi
    else
        fprintf('out simplex')
    end
    
    
    
    %  Simplex two
    f(:,1) = sis_func(d11);
    f(:,2) = sis_func(d21);
    f(:,3) = sis_func(d22);
    T = [(f(:,1)-f(:,3)) (f(:,2)-f(:,3))]
    
    lambda = inv(T)*(-(f(:,3)))
    
    if(lambda >= 0 & sum(lambda)<=1)
        fprintf('no simplexo')
    else
        fprintf('fora simplexo')
    end
    
    %% Defines surface
    %article cazals
    U = [0 0 0 0 0 1 1 1 1 1];
    V = [0 0 0 0 0 1 1 1 1 1];
    k = 4;
    l = 4;
    P        = [0 1/4 2/4 3/4 1;   0  1/4 2/4 3/4 4/4;   0  1/4 2/4 3/4 4/4;   0  1/4 2/4 3/4 4/4;  0 1/4 2/4 3/4 1];
    P(:,:,2) = [0  0   0   0  0;  1/4 1/4 1/4 1/4 1/4;  2/4 2/4 2/4 2/4 2/4;  3/4 3/4 3/4 3/4 3/4;  1  1   1   1  1];
    P(:,:,3) = [0  0   0   0  0;   0   1  -1  -1   0 ;   0  -1   1   1   0 ;   0   1  -1   1   0 ;  0  0   0   0  0];
    P(:,:,4) = [1  1   1   1  1;   1   1   1   1   1 ;   1   1   1   1   1 ;   1   1   1   1   1 ;  1  1   1   1  1];
    %corrects parameters
    P(:,:,1) = transpose(P(:,:,1));
    P(:,:,2) = transpose(P(:,:,2));
    P(:,:,3) = transpose(P(:,:,3));
    P(:,:,4) = transpose(P(:,:,4));
    %homogenous coordinates
    PW(:,:,1)   = P(:,:,4).*P(:,:,1);
    PW(:,:,2)   = P(:,:,4).*P(:,:,2);
    PW(:,:,3)   = P(:,:,4).*P(:,:,3);
    PW(:,:,4)   = P(:,:,4);
    
    %% Surface operations
    %derivatives of the surface (order 1 to 3)
    Du   = derivative_order_cd_NURBS_surface(u,v,U,V,PW,k+1,l+1,2,1);
    Dv   = derivative_order_cd_NURBS_surface(u,v,U,V,PW,k+1,l+1,1,2);
    Duu  = derivative_order_cd_NURBS_surface(u,v,U,V,PW,k+1,l+1,3,1);
    Duv  = derivative_order_cd_NURBS_surface(u,v,U,V,PW,k+1,l+1,2,2);
    Dvv  = derivative_order_cd_NURBS_surface(u,v,U,V,PW,k+1,l+1,1,3);
    Duuu = derivative_order_cd_NURBS_surface(u,v,U,V,PW,k+1,l+1,4,1);
    Duuv = derivative_order_cd_NURBS_surface(u,v,U,V,PW,k+1,l+1,3,2);
    Duvv = derivative_order_cd_NURBS_surface(u,v,U,V,PW,k+1,l+1,2,3);
    Dvvv = derivative_order_cd_NURBS_surface(u,v,U,V,PW,k+1,l+1,1,4);
    Duuuu = derivative_order_cd_NURBS_surface(u,v,U,V,PW,k+1,l+1,5,1);
    Duuuv = derivative_order_cd_NURBS_surface(u,v,U,V,PW,k+1,l+1,4,2);
    Duuvv = derivative_order_cd_NURBS_surface(u,v,U,V,PW,k+1,l+1,3,3);
    Duvvv = derivative_order_cd_NURBS_surface(u,v,U,V,PW,k+1,l+1,2,4);
    Dvvvv = derivative_order_cd_NURBS_surface(u,v,U,V,PW,k+1,l+1,1,5);    
    %first fundamental form
    [E,Eu,Ev,Euu,Euv,Evv,F,Fu,Fv,Fuu,Fuv,Fvv,G,Gu,Gv,Guu,Guv,Gvv] = first_fundamental_form(Du,Dv,Duu,Duv,Dvv,Duuu,Duuv,Duvv,Dvvv);
    %surface normal
    [N,barN,barNu,barNv,barNuu,barNuv,barNvv] = surface_normal(Du,Dv,Duu,Duv,Dvv,Duuu,Duuv,Duvv,Dvvv);
    %second fundamental form
    [e,bare,bareu,barev,bareuu,bareuv,barevv,f,barf,barfu,barfv,barfuu,barfuv,barfvv,g,barg,bargu,bargv,barguu,barguv,bargvv] = second_fundamental_form(Duu,Duv,Dvv,Duuu,Duuv,Duvv,Dvvv,Duuuu,Duuuv,Duuvv,Duvvv,Dvvvv,N,barN,barNu,barNv,barNuu,barNuv,barNvv);
    %auxiliary terms
    [A,Au,Av,Auu,Auv,Avv,B,barB,barBu,barBv,barBuu,barBuv,barBvv,C,barC,barCu,barCv,barCuu,barCuv,barCvv] = auxiliary_ABC(E,Eu,Ev,Euu,Euv,Evv,F,Fu,Fv,Fuu,Fuv,Fvv,G,Gu,Gv,Guu,Guv,Gvv,e,bare,bareu,barev,bareuu,bareuv,barevv,f,barf,barfu,barfv,barfuu,barfuv,barfvv,g,barg,bargu,bargv,barguu,barguv,bargvv);
    [Pu,Puu,Puv,Pv,Pvu,Pvv,Ru,Ruu,Ruv,Rv,Rvu,Rvv,Q,Qu,Qv] = auxiliary_PRQ(A,Au,Av,Auu,Auv,Avv,barB,barBu,barBv,barBuu,barBuv,barBvv,barC,barCu,barCv,barCuu,barCuv,barCvv);
    
    %% System of equations
    y(1) = Q.*(Pu.^2) - Ru.^2;
    y(2) = Q.*(Pv.^2) - Rv.^2;
end