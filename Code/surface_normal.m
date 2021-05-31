function [N,barN,barNu,barNv,barNuu,barNuv,barNvv] = surface_normal(Du,Dv,Duu,Duv,Dvv,Duuu,Duuv,Duvv,Dvvv)
    %cross product
    barN = cross(Du,Dv);                %normal vector
    cross_product_norm  = norm(barN);
    
    %unitary normal vector
    N = barN./cross_product_norm;
    
    %partial derivatives
    barNu  = cross(Duu,Dv) + cross(Du,Duv);
    barNv  = cross(Duv,Dv) + cross(Du,Dvv);
    barNuu = cross(Duuu,Dv) + cross(Duu,Duv).*2 + cross(Du,Duuv);
    barNuv = cross(Duuv,Dv) + cross(Duu,Dvv) + cross(Duv,Duv) + cross(Du,Duvv);
    barNvv = cross(Duvv,Dv) + cross(Duv,Dvv).*2 + cross(Du,Dvvv);
end