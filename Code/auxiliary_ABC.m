function [A,Au,Av,Auu,Auv,Avv,B,barB,barBu,barBv,barBuu,barBuv,barBvv,C,barC,barCu,barCv,barCuu,barCuv,barCvv] = auxiliary_ABC(E,Eu,Ev,Euu,Euv,Evv,F,Fu,Fv,Fuu,Fuv,Fvv,G,Gu,Gv,Guu,Guv,Gvv,e,bare,bareu,barev,bareuu,bareuv,barevv,f,barf,barfu,barfv,barfuu,barfuv,barfvv,g,barg,bargu,bargv,barguu,barguv,bargvv)
    %irrational
    B = 2.*F.*f - G.*e - E.*g;
    C = e.*g - f.^2;
    
    %rational
    A    = E.*G - F.^2;
    barB = 2.*barf.*F - bare.*G - barg.*E;
    barC = bare.*barg - barf.^2;
    
    %partial derivatives 1st order
    Au    = Eu.*G + Gu.*E - 2.*F.*Fu;
    Av    = Ev.*G + Gv.*E - 2.*F.*Fv;
    barBu = 2.*(F.*barfu + Fu.*barf) - (Eu.*barg + E.*bargu) - (Gu.*bare + G.*bareu);
    barBv = 2.*(F.*barfv + Fv.*barf) - (Ev.*barg + E.*bargv) - (Gv.*bare + G.*barev);
    barCu = bareu.*barg + bare.*bargu - 2.*barf.*barfu;
    barCv = barev.*barg + bare.*bargv - 2.*barf.*barfv;
    
    %partial derivatives 2nd order
    Auu    = Euu.*G + 2.*Gu.*Eu + Guu.*E - 2.*(Fu.*Fu + F.*Fuu);
    Auv    = Euv.*G + Eu.*Gv + Guv.*E + Gu.*Ev - 2.*(Fv.*Fu + F.*Fuv);
    Avv    = Evv.*G + 2.*Gv.*Ev + Gvv.*E - 2.*(Fv.*Fv + F.*Fvv);
    barBuu = 2.*(Fu.*barfu + F.*barfuu + Fuu.*barf + Fu.*barfu) - (Euu.*barg + Eu.*bargu + Eu.*bargu + E.*barguu) - (Guu.*bare + Gu.*bareu + Gu.*bareu + G.*bareuu);
    barBuv = 2.*(Fv.*barfu + F.*barfuv + Fuv.*barf + Fu.*barfv) - (Euv.*barg + Eu.*bargv + Ev.*bargu + E.*barguv) - (Guv.*bare + Gu.*barev + Gv.*bareu + G.*bareuv);
    barBvv = 2.*(Fv.*barfv + F.*barfvv + Fvv.*barf + Fv.*barfv) - (Evv.*barg + Ev.*bargv + Ev.*bargv + E.*bargvv) - (Gvv.*bare + Gv.*barev + Gv.*barev + G.*barevv);
    barCuu = bareuu.*barg + bareu.*bargu + bareu.*bargu + bare.*barguu - 2.*barfu.*barfu - 2.*barf.*barfuu;
    barCuv = bareuv.*barg + bareu.*bargv + barev.*bargu + bare.*barguv - 2.*barfv.*barfu - 2.*barf.*barfuv;
    barCvv = barevv.*barg + barev.*bargv + barev.*bargv + bare.*bargvv - 2.*barfv.*barfv - 2.*barf.*barfvv;
end