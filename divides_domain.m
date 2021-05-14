function [testk1, testk2, test2] = divides_domain(umin,umax,vmin,vmax,divide_x,divide_y)
    space_x = (umax - umin)./divide_x;
    space_y = (vmax - vmin)./divide_y;
    space = [space_x space_y];
    
    %starts from [0,0]
    start00 = [0 0];
    start10 = [space_x 0];
    start01 = [0 space_y];
    start11 = [space_x space_y];
    
    countk1  = 0;
    countk2 = 0;
    count2 = 0;
    testk1 = [0 0];
    testk2 = [0 0];
    test2 = [0 0];
    tic
    for i=1:divide_x
       for j=1:divide_y
          square00 = start00 + (i-1).*[space_x 0] + (j-1).*[0 space_y];
          square10 = start10 + (i-1).*[space_x 0] + (j-1).*[0 space_y];
          square01 = start01 + (i-1).*[space_x 0] + (j-1).*[0 space_y];
          square11 = start11 + (i-1).*[space_x 0] + (j-1).*[0 space_y];
          [x,y,z,flag,u] = check_simplexes(square01,square11,square00,square10);
          [x2,y2,z2,flag2,u2] = check_simplexes2(square01,square11,square00,square10);
          if flag2
             [regulafalsi2] = regula_falsi_2d2(x2,y2,z2,3);
             [answer2 n2 status2] = Newton_Sis2(regulafalsi2,1.0e-5,500);
             %[values] = check_critical_condition(answer);
             [values2] = check_umbilic_condition(answer2);
             %if ((values(1)<=1.0e-6 | values(3)<=1.0e-7) & (values(2)<=1.0e-7 | values(4)<=1.0e-7))
             if values2 < 1.0e-6 & values2 > -1.0e-6
                 count2 = count2 + 1;
                 test2(count2,:) = answer2;
             end
          end
          if flag
             [regulafalsi] = regula_falsi_2d(x,y,z,3);
             [answer n status] = Newton_Sis(regulafalsi,1.0e-5,500);
             [values1u,values1v,values2u,values2v] = check_critical_condition(answer);
             if (values1u < 1.0e-5) && (values1u > -1.0e-5) && (values1v < 1.0e-5) && (values1v > -1.0e-5)
                 countk1 = countk1 + 1;
                 testk1(countk1,:) = answer;
             end
             if (values2u < 1.0e-5) && (values2u > -1.0e-5) && (values2v < 1.0e-5) && (values2v > -1.0e-5)
                 countk2 = countk2 + 1;
                 testk2(countk2,:) = answer;
             end
          end
       end
    end
    toc
    hold on
    daspect([1 1 1]);
    scatter(test2(:,1),test2(:,2),'k')
    scatter(testk1(:,1),testk1(:,2),'b')
    scatter(testk2(:,1),testk2(:,2),'r')
    axis([0 1 0 1])
    hold off
end