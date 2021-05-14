function [i] = find_knot_interval(u,U,k)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%Determine knot interval [u_i,u_i+1) u belongs
    %%Input:  Parameter u
    %%        Knot Vector U
    %%        B-Spline degree is k-1
    %%Return: the index i from u_i
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    m = size(U,2); %how many knots?
    if(u == U(m))
        i = m;     %returns last element
    %binary search
    else
       low = k-1;
       high = m;
       i = floor((low+high)/2);
       while(u < U(i+1) || u >= U(i+2))
           if(u < U(i+1))
               high = i;
           else
               low = i;
           end
           i = floor((low+high)/2);
       end
       i = i+1;    %returns index i (for matlab)
    end
end