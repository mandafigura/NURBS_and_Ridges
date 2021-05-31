function [answer] = regula_falsi_2d(the_function,a,b,c,itmax)
    %% Generalized Regula-Falsi Method for 2 dimensions
    %% Description:
    %               finds an approximate value of the zero of a function F
    %               inside a 2-simplex. Algorithm from the article: "Simplicial
    %               Methods for the Solution of Systems of Nonlinear Equations"
    %- Input:
    %         sis_func, the function F
    %         a b c,    each one is a vector of size 2 representing the
    %                   2-simplex coordinates
    %         itmax,    max number of iterations for the method
    %- Output:
    %         answer,   the best approximation after itmax iterations
    %% step (i)
    x(1,:) = a;
    x(2,:) = b;
    x(3,:) = c;
    counter = 0;
    while(counter <= itmax)
        %% step (ii)
        f0 = the_function(x(1,:));
        f1 = the_function(x(2,:));
        f2 = the_function(x(3,:));
        C = [1 1 1; f0' f1' f2'];
        d = [1; 0; 0];
        lambda = lsqnonneg(C,d);
        %% step (iii)
        answer = sum(diag(lambda)*[x(1,:);x(2,:);x(3,:)]);
        fx = the_function(answer);
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