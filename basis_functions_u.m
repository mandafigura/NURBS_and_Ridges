function [N] = basis_functions_u(u,U,n)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%Computes all the basis functions values N_{i,k}(u) til degree n-1.
    %%Based on the recursive definition of the basis functions.
    %%Input:  
    %%        Parameter u
    %%        Knot Vector U
    %%        B-Spline degree is n-1
    %%Return: 
    %%        A matrix N(m,n) with all N_{i,k}(u) values (m is the size of U)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    m = size(U,2);                %how many knots?
    for k = 1:n                   %runs through basis functions degrees
       if k == 1                  %initial case, k = 0
           for i = 1:(m-1)        %runs through knots
               if U(i) <= u & u < U(i+1)
                   N(i,1) = 1;
               else
                   N(i,1) = 0;
               end
           end
       else                       %Recursive definition of N_{i,k}
           for i = 1:(m-k)        %runs through knots
               N(i,k) = 0;  
               if N(i,k-1) ~= 0   %computes equation first term
                   N(i,k) = N(i,k-1).*((u-U(i))./(U(i+k-1)-U(i)));
               end
               if N(i+1,k-1) ~= 0 %computes equation second term
                   N(i,k) = N(i,k) + N(i+1,k-1).*((U(i+k)-u)./(U(i+k)-U(i+1)));
               end
           end
       end
    end
end