clear

%knot vectors of the curve
%U = [0 0 0 0 1 2 3 4 4 4 4];
U = [0 0 0 1 2 3 4 5 6 7 8 9 10 10 10]; %heart
%U = [0 0 0 0 1 2 3 4 4 4 4];
%U = [0 0 0 0 1 2 3 4 5 6 6 6 6];
%U = [0 0 0 0 0 0 1 2 3 5 8 13 21 34 55 89 89 89 89 89 89]
%U = [0 0 0 1 2 3 4 4 5 5 5];
%U = [0 0 0 1 2 2 2];

%degree of the curve
k  = 2;
%derivative order
p = 1;
%point on the curve
point = 0.7

%the domain of the curve
u  = linspace(min(U),max(U),5000);

%Points and Weights that define the curve (heart)
W = [1 1 1 1 1 1 1 1 1 1 1 1];
PX = [2.5 0 0 1 2 2.5 2.5 3 4 5 5 2.5];
PY = [0 1 3 4 3.5 2.7 2.7 3.5 4 3 1 0];

%PX = [0 2 2 3 3 3 4 5 6];
%PY = [1 0 2 2 2 2 0 0 1];

%W = [1 1 1 1 1 1 1 1];
%PX = [1 2 3 4 5 6 7 8];
%PY = [0 1 0 4 1 3 0 2];

%PX = [0 2 2 3 4 5 6 7 8 9 10 11 12 13 14];
%PY = [1 0 2 2 0 0 1 2 0 5 8 2 8 3 -1];

%W =  [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]
%homogenous coordinates
P        = W.*PX;
P(2,:)   = W.*PY;
P(3,:)   = W;

%Curve
    %Homogenous curve
for i = 1:(length(u)-1)
    [SX(i), SY(i), SZ(i)] = parameter_u_b_spline_curve(u(i),U,P,k+1);
end
    %Homogenous projection
for i = 1:(length(u)-1)
    SW(1) = SX(i);
    SW(2) = SY(i);
    SW(3) = SZ(i);
    [H] = H_Perspective(SW);
    SWX(i) = H(1);
    SWY(i) = H(2);
end
%end of curve

%Derivative of the curve
[DSX,DSY] = derivative_order_c_NURBS_curve(point,U,P,k+1,p+1);
derivativeX = 0 + DSX(1);
derivativeY = 0 + DSY(1);
exibition   = 0.5
derivativeX(2) = exibition.*DSX(p+1) + DSX(1);
derivativeY(2) = exibition.*DSY(p+1) + DSY(1);

%start plotting
hold on
grid on;
box on;
%daspect([2 5 1]);
xlabel('eixo x');               %x-axis legend
ylabel('eixo y');               %y-axis legend
txt1 = ['Pontos de Controle'];
%set(gca,'visible','off')
txt2 = ['Curva Spline'];
txt0 = ['Derivada'];
%plot(SWX, SWY,'r-', 'LineWidth',1.7, 'DisplayName',txt2);
plot(SWX, SWY,'r-', 'LineWidth',2.7, 'DisplayName',txt2);               %curve
plot(PX,PY,'b-*','LineWidth',1.2,'DisplayName',txt1);                   %Points
plot(derivativeX,derivativeY,'k-o','LineWidth',1.2,'DisplayName',txt0); %Derivative
%txt3 = ['Pontos de Ruptura'];
%scatter(KX(k+2:length(U)-k-1), KY(k+2:length(U)-k-1),'k', 'DisplayName',txt3);
%legend('Location','eastoutside');
%legend('Location','northeast');
%axis([-0.2 5.2 -0.2 4.7])
hold off