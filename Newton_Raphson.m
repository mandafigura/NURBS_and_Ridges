function [x n status] = Newton_Raphson(the_function,the_function_jacobian,x,tol,nmax)
   %% Newton-Raphson method for solving equations
   %- Input:
   %         sis_func,  the function F
   %         Dsis_func, the jacobian JF
   %         x,         initial guess
   %         tol,       tolerance error
   %         nmax,      max number of iterations
   %- Output:
   %         x,         the solution
   %         n,         number of iterations
   %         status,    in case of convergence: true, didn't converge: false
   
   %% sets beggining
   n = 0;
   error = 2*tol*norm(x);
   
   %% Method
   while n < nmax && error >= tol*norm(x)
       f  = the_function(x);
       df = the_function_jacobian(x);
       y  = -df\f';
       error = norm(y);
       x    = x+y';
       n    = n+1;
   end
   
   %% Status
   if n < nmax
       status = true;
   else
       status = false;
   end
end