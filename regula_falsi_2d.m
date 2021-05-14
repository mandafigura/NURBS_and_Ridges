function [answer] = regula_falsi_2d(a,b,c,itmax)
    %% step (i)
    x(1,:) = a;
    x(2,:) = b;
    x(3,:) = c;
    counter = 0;
    while(counter <= itmax)
        %% step (ii)
        f0 = sis_func(x(1,:));
        f1 = sis_func(x(2,:));
        f2 = sis_func(x(3,:));
        C = [1 1 1; f0' f1' f2'];
        d = [1; 0; 0];
        lambda = lsqnonneg(C,d);
        %% step (iii)
        answer = sum(diag(lambda)*[x(1,:);x(2,:);x(3,:)]);
        fx = sis_func(answer);
        if(fx == [0 0])
            return
        end
        %% step (iv)
        [xi] = baricentric_coordinates_2d(f0,f1,f2,fx);
        %% step (v)
        test = lambda./xi;
        for i=1:3
            [m,index] = min(test);
            if xi(index) > 0
                break;
            else
                test(index) = Inf;
            end
        end
        %% step (vi)
        x(index,:) = answer;
        counter = counter + 1;
    end
end