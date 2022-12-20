function [ridge_points] = predictor_corrector(x,h,hmin,hmax,bark,maxk,bard,bara,tol,maxit,Func,DFunc,Tangent,maxtrials,umin,umax,vmin,vmax,criticalist,umbiliclist,index)
    %hmin = 1.0e-15;
    %maxtrials = 15;
    %% Description:;
    %               Algorithm based on the book " Introduction to Numerical
    %               Cotinuation Methods" (Eugene L. Allgower and Kurt Georg), 1990.
    %               Algorithm (6.1.10), page 48.
    % Input:
    %         x,     initial point such that Func(x) = 0 (row vector)
    %         h,     initial step lenght                 (scalar)
    %         hmin,  minimum steplenght                  (scalar)
    %         hmax,  maximum steplenght                  (scalar)
    %         bark,  nominal contraction rate            (scalar)
    %         mark,  maximum tolerated contraction rate  (scalar)
    %         bard,  nominal distance to curve           (scalar)
    %         bara,  nominal angle                       (scalar)
    %         tol,   error colerance for Newton Method on Corrector function
    %         maxit, max number of irerations before corrector method convergence
    %         Func,  the function to be traced
    %         Dfunc, derivative of the function to be traced
    %         maxtrials, maximum number of trials for adjusting h... forces a step with hmin size
    %         umin,umax,vmin,vmax,       where the B-Splines are defined
    %         criticalist,umbiliclist,   seedpoints
    steps = 0;
    ridge_points = x;
    force_next_step = 0;
    keepgoing = true;
    frustration = 0;
    prev_tangent = [Inf Inf];
    while keepgoing% && steps < 400
        DHx  = DFunc(x);            % evaluates derivative DH at x
        tDHx = Tangent(DHx);  % gets tangent of function H at x
        if prev_tangent ~= [Inf Inf]        % tangentes estão na mesma direção?
           
           %prev_tangent
           %tDHx
           direction = dot(tDHx,prev_tangent);
           if direction < 0                 % verifica se formam um angulo agudo
               %DHx
               %tDHx
               tDHx = (-1)*tDHx;
               %steps
               fprintf('K')
           end
        end
        %% Adjusts tangent direction
        %test_matrix = [DHx;tDHx];
        %if det(test_matrix)<0; tDHx = (-1)*tDHx; end;
        %% Predictor Step
        y = x + h*tDHx;             % a row vector
        [keepgoing,flag] = predictor_stop_condition(y,umin,umax,vmin,vmax,criticalist,umbiliclist,hmax,index);
        if ~keepgoing; return; end; % should we stop the trace?
        %% Corrector loop
        [z,converge,delta,kontr,alpha] = corrector(y,tDHx,maxk,maxit,tol,Func,DFunc,Tangent);
        if converge
            fprintf('A');
            %% Steplength adaptation
            f = max([sqrt(kontr/bark) sqrt(delta/bard) (alpha/bara)]);
            f = max([min([f 2]) (1/2)]); % deceleration factor   
            h = min([max([hmin h/f]) hmax]); % steplength adaptation
            if f<2                       % new point is accepted
                x = z;
                steps = steps + 1;
                ridge_points((steps+1),:) = x;
                fprintf('\n')
                force_next_step = 0;     % starts a new step
                frustration = 0;
                prev_tangent = tDHx;
            else
                force_next_step = force_next_step + 1; % counts how many failed step sizes
                if force_next_step > maxtrials    %forces a new step
                    x = z;
                    steps = steps + 1;
                    ridge_points((steps+1),:) = x;
                    fprintf('M\n')
                    %force_next_step = 0;
                    frustration = 0;
                    h = hmin;
                    prev_tangent = tDHx;
                end
            end
        else        %if the corrector doen't converges properly, adjusts step size
            h = min([max([hmin h/2]) hmax]);
            frustration = frustration + 1;
            fprintf('F')
            if frustration > 20; keepgoing = false; end;
        end
    end
end