function [Pu,Puu,Puv,Pv,Pvu,Pvv,Ru,Ruu,Ruv,Rv,Rvu,Rvv,Q,Qu,Qv] = auxiliary_PRQ(A,Au,Av,Auu,Auv,Avv,barB,barBu,barBv,barBuu,barBuv,barBvv,barC,barCu,barCv,barCuu,barCuv,barCvv)
    Pu = ((-A.^(-1.5)).*barBu  +  (1.5).*(A.^(-2.5)).*Au.*barB)./2;
    %partial derivatives of Pu
    Puu = (3/4).*((A.^(-2.5)).*Au.*barBu) - (1/2).*((A.^(-1.5)).*barBuu) - (15/8).*((A.^(-3.5)).*Au.*Au.*barB) + (3/4).*((A.^(-2.5)).*Auu.*barB) + (3/4).*((A.^(-2.5)).*Au.*barBu);
    Puv = (3/4).*((A.^(-2.5)).*Av.*barBu) - (1/2).*((A.^(-1.5)).*barBuv) - (15/8).*((A.^(-3.5)).*Au.*Av.*barB) + (3/4).*((A.^(-2.5)).*Auv.*barB) + (3/4).*((A.^(-2.5)).*Au.*barBv);
    
    Pv = ((-A.^(-1.5)).*barBv  +  (1.5).*(A.^(-2.5)).*Av.*barB)./2;
    %partial derivatives of Pv
    Pvu = (3/4).*((A.^(-2.5)).*Au.*barBv) - (1/2).*((A.^(-1.5)).*barBuv) - (15/8).*((A.^(-3.5)).*Au.*Av.*barB) + (3/4).*((A.^(-2.5)).*Auv.*barB) + (3/4).*((A.^(-2.5)).*Av.*barBu);
    Pvv = (3/4).*((A.^(-2.5)).*Av.*barBv) - (1/2).*((A.^(-1.5)).*barBvv) - (15/8).*((A.^(-3.5)).*Av.*Av.*barB) + (3/4).*((A.^(-2.5)).*Avv.*barB) + (3/4).*((A.^(-2.5)).*Av.*barBv);

    
    Ru = ((A.^(-1.5)).*barBu.*barB  -  2.*(A.^(-0.5)).*barCu  +  4.*(A^(-1.5)).*Au.*barC  -  (1.5).*(A.^(-2.5)).*Au.*(barB.^2))./2;
    %partial derivatives of Ru
    Ruu = (-3/4).*((A.^(-2.5)).*Au.*barB.*barBu) + (1/2).*((A.^(-1.5)).*barB.*barBuu) + (1/2).*((A.^(-1.5)).*barBu.*barBu) + (1/2).*((A.^(-1.5)).*Au.*barCu) - ((A.^(-0.5)).*barCuu) - (3.*(A.^(-2.5)).*Au.*Au.*barC) + (2.*(A.^(-1.5)).*Auu.*barC) + (2.*(A.^(-1.5)).*Au.*barCu) + (15/8).*((A.^(-3.5)).*Au.*Au.*barB.*barB) - (3/4).*((A.^(-2.5)).*Auu.*barB.*barB) - (3/2).*((A.^(-2.5)).*Au.*barB.*barBu);
    Ruv = (-3/4).*((A.^(-2.5)).*Av.*barB.*barBu) + (1/2).*((A.^(-1.5)).*barB.*barBuv) + (1/2).*((A.^(-1.5)).*barBu.*barBv) + (1/2).*((A.^(-1.5)).*Av.*barCu) - ((A.^(-0.5)).*barCuv) - (3.*(A.^(-2.5)).*Au.*Av.*barC) + (2.*(A.^(-1.5)).*Auv.*barC) + (2.*(A.^(-1.5)).*Au.*barCv) + (15/8).*((A.^(-3.5)).*Au.*Av.*barB.*barB) - (3/4).*((A.^(-2.5)).*Auv.*barB.*barB) - (3/2).*((A.^(-2.5)).*Au.*barB.*barBv);
    
    Rv = ((A.^(-1.5)).*barBv.*barB  -  2.*(A.^(-0.5)).*barCv  +  4.*(A^(-1.5)).*Av.*barC  -  (1.5).*(A.^(-2.5)).*Av.*(barB.^2))./2;
    %partial derivatives of Rv
    Rvu = (-3/4).*((A.^(-2.5)).*Au.*barB.*barBv) + (1/2).*((A.^(-1.5)).*barB.*barBuv) + (1/2).*((A.^(-1.5)).*barBu.*barBv) + (1/2).*((A.^(-1.5)).*Au.*barCv) - ((A.^(-0.5)).*barCuv) - (3.*(A.^(-2.5)).*Au.*Av.*barC) + (2.*(A.^(-1.5)).*Auv.*barC) + (2.*(A.^(-1.5)).*Av.*barCu) + (15/8).*((A.^(-3.5)).*Au.*Av.*barB.*barB) - (3/4).*((A.^(-2.5)).*Auv.*barB.*barB) - (3/2).*((A.^(-2.5)).*Av.*barB.*barBu);
    Rvv = (-3/4).*((A.^(-2.5)).*Av.*barB.*barBv) + (1/2).*((A.^(-1.5)).*barB.*barBvv) + (1/2).*((A.^(-1.5)).*barBv.*barBv) + (1/2).*((A.^(-1.5)).*Av.*barCv) - ((A.^(-0.5)).*barCvv) - (3.*(A.^(-2.5)).*Av.*Av.*barC) + (2.*(A.^(-1.5)).*Avv.*barC) + (2.*(A.^(-1.5)).*Av.*barCv) + (15/8).*((A.^(-3.5)).*Av.*Av.*barB.*barB) - (3/4).*((A.^(-2.5)).*Avv.*barB.*barB) - (3/2).*((A.^(-2.5)).*Av.*barB.*barBv);
    
    Q  = (barB.^2) - 4.*A.*barC;
    %partial derivatives of Q
    Qu = 2.*barB.*barBu - 4.*(Au.*barC + A.*barCu);
    Qv = 2.*barB.*barBv - 4.*(Av.*barC + A.*barCv);
end