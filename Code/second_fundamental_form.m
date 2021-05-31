function [e,bare,bareu,barev,bareuu,bareuv,barevv,f,barf,barfu,barfv,barfuu,barfuv,barfvv,g,barg,bargu,bargv,barguu,barguv,bargvv] = second_fundamental_form(Duu,Duv,Dvv,Duuu,Duuv,Duvv,Dvvv,Duuuu,Duuuv,Duuvv,Duvvv,Dvvvv,N,barN,barNu,barNv,barNuu,barNuv,barNvv)
    %second fundamental form
    e = dot(Duu,N);
    f = dot(Duv,N);
    g = dot(Dvv,N);
    
    %rational second fundamental form
    bare  = dot(Duu,barN);
    bareu = dot(barNu,Duu) + dot(barN,Duuu);    %partial derivative u
    barev = dot(barNv,Duu) + dot(barN,Duuv);    %partial derivative v
    barf  = dot(Duv,barN);
    barfu = dot(barNu,Duv) + dot(barN,Duuv);    %partial derivative u
    barfv = dot(barNv,Duv) + dot(barN,Duvv);    %partial derivative v
    barg  = dot(Dvv,barN);
    bargu = dot(barNu,Dvv) + dot(barN,Duvv);    %partial derivative u
    bargv = dot(barNv,Dvv) + dot(barN,Dvvv);    %partial derivative v
    
    %second order partial derivatives
    bareuu = dot(barNuu,Duu) + dot(barNu,Duuu) + dot(barNu,Duuu) + dot(barN,Duuuu);
    bareuv = dot(barNuv,Duu) + dot(barNu,Duuv) + dot(barNv,Duuu) + dot(barN,Duuuv);
    barevv = dot(barNvv,Duu) + dot(barNv,Duuv) + dot(barNv,Duuv) + dot(barN,Duuvv);
    barfuu = dot(barNuu,Duv) + dot(barNu,Duuv) + dot(barNu,Duuv) + dot(barN,Duuuv);
    barfuv = dot(barNuv,Duv) + dot(barNu,Duvv) + dot(barNv,Duuv) + dot(barN,Duuvv);
    barfvv = dot(barNvv,Duv) + dot(barNv,Duvv) + dot(barNv,Duvv) + dot(barN,Duvvv);
    barguu = dot(barNuu,Dvv) + dot(barNu,Duvv) + dot(barNu,Duvv) + dot(barN,Duuvv);
    barguv = dot(barNuv,Dvv) + dot(barNu,Dvvv) + dot(barNv,Duvv) + dot(barN,Duvvv);
    bargvv = dot(barNvv,Dvv) + dot(barNv,Dvvv) + dot(barNv,Dvvv) + dot(barN,Dvvvv);
end