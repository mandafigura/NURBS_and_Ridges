function [DN] = diff_basis_func(n,u,U,c,i)
    %% ---BASED ON THE NURBS BOOK ALGORITHM 2.5---
    %% Computes the (c-1)th derivative of basis functions of degree n-1
    % Input:  
    %         Parameter u
    %         Knot Vector U
    %         B-Spline degree is n-1
    %         Diff order is c-1
    %         N_i-1,n(u) to be computed (i)
    % Return: 
    %         Derivative array DN

    %% Computes basis functions
    %local property
    if(u < U(i) | u >= U(i+n))
        for k=1:c
            DN(k) = 0;
        end
    else
        %zero-th degree of basis func
        for j=0:(n-1)
            if(u >= U(i+j) & u < U(i+j+1))
                N(j+1,1) = 1;
            else
                N(j+1,1) = 0;
            end
        end
        %computes others degrees of basis func
        for k=2:n
            if(N(1,k-1) == 0)
                saved = 0;
            else
                saved = ((u-U(i)).*N(1,k-1))./(U(i+k-1)-U(i));
            end
            for j=0:(n-k)
               Uleft = U(i+j+1);
               Uright = U(i+j+k);
               if(N(j+2,k-1) == 0)
                   N(j+1,k) = saved;
                   saved = 0;
               else
                   temp = N(j+2,k-1)./(Uright-Uleft);
                   N(j+1,k) = saved+((Uright - u).*temp);
                   saved = (u-Uleft).*temp;
               end
            end
        end
        %% Computes derivatives
        %degree 0 derivative
        DN(1) = N(1,n);
        %compute derivatives of each degree til degree c-1
        for k=2:c
            %load appropriate column
            for j=1:k
                ND(j) = N(j,n-k+1);
            end
            for jj=1:(k-1)
               if(ND(1) == 0)
                   saved = 0;
               else
                   saved = ND(1)./(U(i+n-k+jj)-U(i));
               end
               for j=0:(k-jj-1)
                  Uleft = U(i+j+1);
               %%%%alteracao para evitar right nao definido
                  test = i+j+n+jj-1;
                  if(test <= length(U))
                    Uright = U(i+j+n+jj-1);  %%%linha original
                  end
               %%%%%%%%fim alteração
                  if(ND(j+2) == 0)
                      ND(j+1) = (n-k+jj).*saved;
                      saved  = 0;
                  else
                  %%%%%alteração para nao dar inf (divisao por 0)
                      if((Uright-Uleft)~=0)
                          temp = ND(j+2)./(Uright-Uleft);  %%linha original
                      end
                  %%%%%fim alteracao
                     ND(j+1) = (n-k+jj).*(saved-temp);
                     saved = temp;
                  end
               end
            end
            DN(k) = ND(1);
        end    
    end
end