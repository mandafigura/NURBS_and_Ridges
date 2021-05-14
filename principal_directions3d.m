function [T1,T2] = principal_directions3d(u,v)
    %return normalized principal directions vectors
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
    %principal directions 2d
    [t1,t2] = principal_directions2d(u,v);
    
    T1 = t1(1,1).*Du + t1(2,1).*Dv;
    T2 = t2(1,1).*Du + t2(2,1).*Dv;
    
    T1 = T1./norm(T1);
    T2 = T2./norm(T2);
end