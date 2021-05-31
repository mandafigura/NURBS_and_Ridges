function [criticalk1, criticalk2, umbilics] = finds_seed_points(umin,umax,vmin,vmax,subdivision_x,subdivision_y)
    warning off;
    spacing_x = (umax - umin)./subdivision_x;
    spacing_y = (vmax - vmin)./subdivision_y;
    
    %% Starts dividing the domain from point [0,0]
    %  Makes rectangles in the following order:
    %
    %        3 ----- 4
    %        |       |
    %        |       |
    %        1 ----- 2
    %
    %% The first retangle in the domain
    start00 = [0 0];
    start10 = [spacing_x 0];
    start01 = [0 spacing_y];
    start11 = [spacing_x spacing_y];
    
    %% Initial auxiliary values
    count_k1 = 0;          % How many k1 critical points?
    count_k2 = 0;          % How many k2 critical points?
    count_umbilics = 0;    % How many umbilics?
    %% Starts variables
    criticalk1 = [0 0];
    criticalk2 = [0 0];
    umbilics   = [0 0];
    %% Values
    tolerance_newton   = 1.0e-5;
    tolerance_critial  = 1.0e-5;
    tolerance_umbilics = 1.0e-6;
    iterations_regulafalsi = 5;
    iterations_newton = 100;
    
    tic
    for i=1:subdivision_x
       for j=1:subdivision_y
          square00 = start00 + (i-1).*[spacing_x 0] + (j-1).*[0 spacing_y];
          square10 = start10 + (i-1).*[spacing_x 0] + (j-1).*[0 spacing_y];
          square01 = start01 + (i-1).*[spacing_x 0] + (j-1).*[0 spacing_y];
          square11 = start11 + (i-1).*[spacing_x 0] + (j-1).*[0 spacing_y];
          [cp_x,cp_y,cp_z,flag_critical_points] = check_simplexes(@critical_points,square01,square11,square00,square10);
          [up_x,up_y,up_z,flag_umbilics]        = check_simplexes(@umbilic_points,square01,square11,square00,square10);
          if flag_umbilics
             [regulafalsi_umbilics] = regula_falsi_2d(@umbilic_points,up_x,up_y,up_z,iterations_regulafalsi);
             [is_it_umbilic n_umbilic status_umbilic] = Newton_Raphson(@umbilic_points,@umbilic_points_jacobian,regulafalsi_umbilics,tolerance_newton,iterations_newton);
             %[values] = check_critical_condition(answer);
             [discriminant] = check_umbilic_condition(is_it_umbilic);
             %if ((values(1)<=1.0e-6 | values(3)<=1.0e-7) & (values(2)<=1.0e-7 | values(4)<=1.0e-7))
             if (discriminant < tolerance_umbilics) & (discriminant > ((-1).*tolerance_umbilics))
                 count_umbilics = count_umbilics + 1;
                 umbilics(count_umbilics,:) = is_it_umbilic;
             end
          end
          if flag_critical_points
             [regulafalsi_critial_points] = regula_falsi_2d(@critical_points,cp_x,cp_y,cp_z,iterations_regulafalsi);
             [is_it_critical n_critical status_critical] = Newton_Raphson(@critical_points,@critical_points_jacobian,regulafalsi_critial_points,tolerance_newton,iterations_newton);
             [kappa_1u,kappa_1v,kappa_2u,kappa_2v] = check_critical_condition(is_it_critical);
             % is it a blue point (maximal, kappa_1)?
             if (kappa_1u < tolerance_critial) && (kappa_1u > -tolerance_critial) && (kappa_1v < tolerance_critial) && (kappa_1v > -tolerance_critial)
                 count_k1 = count_k1 + 1;
                 criticalk1(count_k1,:) = is_it_critical;
             end
             % is it a red point (minimal, kappa_2)?
             if (kappa_2u < tolerance_critial) && (kappa_2u > -tolerance_critial) && (kappa_2v < tolerance_critial) && (kappa_2v > -tolerance_critial)
                 count_k2 = count_k2 + 1;
                 criticalk2(count_k2,:) = is_it_critical;
             end
          end
       end
    end
    toc
    hold on
    daspect([1 1 1]);
    scatter(umbilics(:,1),umbilics(:,2),'k')
    scatter(criticalk1(:,1),criticalk1(:,2),'b')
    scatter(criticalk2(:,1),criticalk2(:,2),'r')
    axis([0 1 0 1])
    hold off
end