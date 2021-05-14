clear                          %clean

%%%%%%%---KNOT VECTOR
%U = [0 0 1 3 3 4 5];           
%U = [0 1 2 3 4 5];
%U = [ 0 0 0 0 1 1 1 1];
%U = [0 1 2 3 4 5 6 7 8 9];
%U = [0 0 0 0 2 4 6 8 8 8 8];
%U = [0 0 0 0 1 2 3 4 4 4 4];
%U = [1 2 3 4 5];
U = [0 0 0 1 2 2 2];
%U = [0 0 0 1 2 3 4 4 5 5 5];
%U = [0 0 0 0 0 0 1 2 3 5 8 13 21 34 55 89 89 89 89 89 89]

%%%%%%%---DEGREE OF B-SPLINE
k = 2; 

%%%%%%%---GENERATE DOMAIN INTERVAL
%u = linspace(min(U)+k+1,max(U)-k,500);
u = linspace(min(U),max(U),500);

%%%%%%%---DERIVATIVE ORDER
g = 1;

for j = 1:(length(u)-1)            %%%RUNS THROUGH THE DOMAIN
    %computes the derivative of all degrees of N_i,k
    for i=1:(length(U)-k-1)
       B(i,:) =  diff_basis_func_on_u(k+1,u(j),U,g+1,i);
    end
    D(:,j) = B(:,g+1); %derivative of order g for each N_i,k(u)
end

%plots basis functions derivative
figure(g+1)
hold on;
box  on;
grid on;
%daspect([1 2 1]);            %proportion
xlabel('eixo x');             %legend x axis
ylabel('eixo y');             %legend y axis

%prints all the basis functions derivatives of order g
for i = 1:(length(U)-k-1)         
    %legend for each function
    txt = ['N^{(',num2str(g) ,')}_{',num2str(i-1), ',', num2str(k), '}'];
    %plots
    plot(u(1:length(u)-1),D(i,:),'DisplayName',txt,'LineWidth',1.7);
end
txt1 = ['Derivatives of order ', num2str(g),' for basis functions of degree ', num2str(k)];
title(txt1);
legend('Location','eastoutside'); %where to show legend
%legend('Location','best');
%legend('Location','northeast');
hold off;