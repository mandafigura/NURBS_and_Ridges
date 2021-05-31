function [gradient] = ridge_gradient(x,which_ridge)
    point = the_coeficients_of_a_point_on_surface(x);
    
    if(strcmp(which_ridge,'blue'))
        k1   = (-(point.B) + sqrt((point.B).^2 - 4.*(point.A).*(point.C)))./(2.*(point.A));
        k1u  = (point.Pu) + ((point.Ru)./sqrt(point.Q));
        k1v  = (point.Pv) + ((point.Rv)./sqrt(point.Q));
        k1uu = (point.Puu) + ((((point.Ruu).*sqrt((point.Q))) - ((point.Ru).*((point.Qu)./(2.*sqrt((point.Q))))))./(point.Q));
        k1uv = (point.Puv) + ((((point.Ruv).*sqrt((point.Q))) - ((point.Ru).*((point.Qv)./(2.*sqrt((point.Q))))))./(point.Q));
        k1vu = (point.Pvu) + ((((point.Rvu).*sqrt((point.Q))) - ((point.Rv).*((point.Qu)./(2.*sqrt((point.Q))))))./(point.Q));
        k1vv = (point.Pvv) + ((((point.Rvv).*sqrt((point.Q))) - ((point.Rv).*((point.Qv)./(2.*sqrt((point.Q))))))./(point.Q));
        t1x  = k1.*(point.F) - (point.f);
        t1y  = (point.e) - k1.*(point.E);
        t1xu = -((((point.barfu).*sqrt(point.A)) - ((point.barf).*((point.Au)./(2.*sqrt(point.A)))))./(point.A)) + (k1u.*(point.F)) + (k1.*(point.Fu));
        t1yu =  ((((point.bareu).*sqrt(point.A)) - ((point.bare).*((point.Au)./(2.*sqrt(point.A)))))./(point.A)) - (k1u.*(point.E)) - (k1.*(point.Eu));
        t1xv = -((((point.barfv).*sqrt(point.A)) - ((point.barf).*((point.Av)./(2.*sqrt(point.A)))))./(point.A)) + (k1v.*(point.F)) + (k1.*(point.Fv));
        t1yv =  ((((point.barev).*sqrt(point.A)) - ((point.bare).*((point.Av)./(2.*sqrt(point.A)))))./(point.A)) - (k1v.*(point.E)) - (k1.*(point.Ev));
        grad1u = (k1uu.*t1x) + (k1u.*t1xu) + (k1vu.*t1y) + (k1v.*t1yu);
        grad1v = (k1uv.*t1x) + (k1u.*t1xv) + (k1vv.*t1y) + (k1v.*t1yv);
        gradient = [grad1u grad1v];
    elseif(strcmp(which_ridge,'red'))
        k2   = (-(point.B) - sqrt((point.B).^2 - 4.*(point.A).*(point.C)))./(2.*(point.A));
        k2u  = (point.Pu) - ((point.Ru)./sqrt(point.Q));
        k2v  = (point.Pv) - ((point.Rv)./sqrt(point.Q));
        k2uu = (point.Puu) - ((((point.Ruu).*sqrt((point.Q))) - ((point.Ru).*((point.Qu)./(2.*sqrt((point.Q))))))./(point.Q));
        k2uv = (point.Puv) - ((((point.Ruv).*sqrt((point.Q))) - ((point.Ru).*((point.Qv)./(2.*sqrt((point.Q))))))./(point.Q));
        k2vu = (point.Pvu) - ((((point.Rvu).*sqrt((point.Q))) - ((point.Rv).*((point.Qu)./(2.*sqrt((point.Q))))))./(point.Q));
        k2vv = (point.Pvv) - ((((point.Rvv).*sqrt((point.Q))) - ((point.Rv).*((point.Qv)./(2.*sqrt((point.Q))))))./(point.Q));
        t2x  = k2.*(point.F) - (point.f);
        t2y  = (point.e) - k2.*(point.E);
        t2xu = -((((point.barfu).*sqrt(point.A)) - ((point.barf).*((point.Au)./(2.*sqrt(point.A)))))./(point.A)) + (k2u.*(point.F)) + (k2.*(point.Fu));
        t2yu =  ((((point.bareu).*sqrt(point.A)) - ((point.bare).*((point.Au)./(2.*sqrt(point.A)))))./(point.A)) - (k2u.*(point.E)) - (k2.*(point.Eu));
        t2xv = -((((point.barfv).*sqrt(point.A)) - ((point.barf).*((point.Av)./(2.*sqrt(point.A)))))./(point.A)) + (k2v.*(point.F)) + (k2.*(point.Fv));
        t2yv =  ((((point.barev).*sqrt(point.A)) - ((point.bare).*((point.Av)./(2.*sqrt(point.A)))))./(point.A)) - (k2v.*(point.E)) - (k2.*(point.Ev));
        grad2u = (k2uu.*t2x) + (k2u.*t2xu) + (k2vu.*t2y) + (k2v.*t2yu);
        grad2v = (k2uv.*t2x) + (k2u.*t2xv) + (k2vv.*t2y) + (k2v.*t2yv);
        gradient = [grad2u grad2v];
    end

end