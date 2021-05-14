function [x2 n2 status2] = Newton_Sis2(x2,tol,nmax)
   n2 = 0;
   erro = 2*tol*norm(x2);
   
   while n2 < nmax && erro >= tol*norm(x2)
       f  = sis_func2(x2);
       df = Dsis_func2(x2);
       y  = -df\f';
       erro = norm(y);
       x2    = x2+y';
       n2    = n2+1;
   end
   
   if n2 < nmax
       status2 = 1;
   else
       status2 = -1;
   end
   
   return

end

