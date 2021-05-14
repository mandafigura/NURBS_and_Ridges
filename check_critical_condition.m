function [values1u,values1v,values2u,values2v] = check_critical_condition(x)
    %% Algorithm description
    %--Generates curvature critical points system of equations
    u = x(1);
    v = x(2);
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
    values1u = Pu + (Ru./sqrt(Q));
    values1v = Pv + (Rv./sqrt(Q));
    values2u = Pu - (Ru./sqrt(Q));
    values2v = Pv - (Rv./sqrt(Q));
end