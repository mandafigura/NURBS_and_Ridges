clear
%Super Simples
%U = [0 0 0 1 1 1]
%superfície com 2 picos
U  = [0 0 0 1 2 3 3 3];
%superfície simples
%U  = [0 0 0 1 2 3 3 3]
%superfície com dobra
%U   = [0 0 0 1 1 1]
u  = linspace(min(U),max(U),100);
%Super Simples
%V = [0 0 0 0 1 1 1 1]
%superfície com 2 picos
V  = [1 1 1 1 2 3 4 5 5 5 5];
%superfície simples
%V  = [0 0 0 0 1 1 1 1]
%superfície com dobra
%V   = [0 0 0 0 1 2 3 4 4 4 4]
v  = linspace(min(V),max(V),100);
k  = 2;
l  = 3;
p = 0;
q = 1;
point = [0.6 4.5];

%Super Simples
%P         = [0 1 2 6; 0 1 2 6; 0 1 2 6];
%P(:,:,2)  = [0 0 0 0; 1 1 1 1; 2 2 2 2];
%P(:,:,3)  = [0 2 4 0; 0 2 4 0; 0 2 4 0];
%superfície com 2 picos
P        = [1 2 3 4 5 6 7; 1 2 3 4 5 6 7; 1 2 3 4 5 6 7; 1 2 3 4 5 6 7; 1 2 3 4 5 6 7];
P(:,:,2) = [0 0 0 0 0 0 0; 1 1 1 1 1 1 1; 2 2 2 2 2 2 2; 3 3 3 3 3 3 3; 4 4 4 4 4 4 4];
P(:,:,3) = [4 4 4 4 4 4 4; 4 4 4 4 4 9 4; 4 4 4 0 4 4 4; 4 9 4 4 4 4 4; 4 4 4 4 4 4 4];

%superfície simples
%P         = [1 2 3 4; 1 2 3 4; 1 2 3 4; 1 2 3 4; 1 2 3 4]
%P(:,:,2)  = [0 0 0 0; 1 1 1 1; 2 2 2 2; 3 3 3 3; 4 4 4 4]
%P(:,:,3)  = [0 0 0 0; 1 1 1 1; 4 4 4 4; 1 1 1 1; 0 0 0 0]

%QQQ = [1 4 0; 1 0 0; 4 0 0; 4 4 0; 4 2 4; 1 2 4]
%[kk,av] = convhulln(QQQ);
%trisurf(kk,QQQ(:,1),QQQ(:,2), QQQ(:,3),'EdgeColor', 'none', 'FaceColor','magenta', 'FaceAlpha', .1)


%superfície com dobra
%P          = [0 1 3 3 3 5 5; 0 1 3 3 3 5 5; 0 1 3 3 3 5 5]
%P(:,:,2)   = [0 0 0 0 0 0 0; 3 3 3 3 3 3 3; 6 6 6 6 6 6 6]
%P(:,:,3)   = [2 4 0 0 0 0 4; 2 4 0 0 0 0 4; 2 4 0 0 0 0 4]


%% Superfície
for i = 1:(length(u)-1)
    for j = 1:(length(v)-1)
        S  = surface_uv_points(k+1,l+1,P,u(i),v(j),U,V);
        SX(i,j) = S(1);
        SY(i,j) = S(2);
        SZ(i,j) = S(3);
    end
end


%% Derivada
DDS = surface_uv_points(k+1,l+1,P,point(1),point(2),U,V)
DD  = derivative_order_cd_bspline_surface(point(1),point(2),U,V,P,k+1,l+1,p+1,q+1)

derivativeX(1) = 0 + DDS(1);
derivativeY(1) = 0 + DDS(2);
derivativeZ(1) = 0 + DDS(3);
exibition   = 1
derivativeX(2) = exibition.*DD(1) + DDS(1)
derivativeY(2) = exibition.*DD(2) + DDS(2)
derivativeZ(2) = exibition.*DD(3) + DDS(3)

%RRR = [0 -1 0; 1 0 0; 0 0 1]
%teste = RRR*transpose(DD)
%derivativeX(2) = teste(1)
%derivativeY(2) = teste(2)
%derivativeZ(2) = teste(3)


%for i = 1:(length(U)-1)
%   [KX(i), KY(i)] = Ponto_u_da_Curva(k+1,PX,PY,U(i),U);
%end

hold on
%a = [-2 -4 0.8];
%a = [-2 -4 0.8]
%a = [-1.8 -4 1]
a = [-4 -1.8 1]
[caz,cel] = view(a)
grid on;
%box on;
%daspect([2 5 1]);
xlabel('eixo x');               % legenda no eixo horizontal
ylabel('eixo y');               % legenda no eixo vertical
zlabel('eixo z');               % legenda no eixo vertical
txt1 = ['Pontos de Controle'];
%set(gca,'visible','off')
txt2 = ['Superfície NUBS'];

%for i = 1:(length(u)-1)
 %   plot3(SX(i,:), SY(i,:), SZ(i,:),'Color',[i/length(u) 0.51 .33], 'LineWidth',1.7, 'DisplayName',txt2);
  %  %plot3(SX(i,:), SY(i,:), SZ(i,:),'Color',[i/length(u) i/length(u) i/length(u)], 'LineWidth',1.7, 'DisplayName',txt2);
  %  %plot3(SX(i,:), SY(i,:), SZ(i,:),'r', 'LineWidth',1.7, 'DisplayName',txt2);
%end
%for j = 1:(length(v)-1)
 %   g(j) = plot3(SX(:,j), SY(:,j), SZ(:,j),'Color',[.77 .22  j/length(v)], 'LineWidth',1.7, 'DisplayName',txt2);
 %   %g(j) = plot3(SX(:,j), SY(:,j), SZ(:,j),'Color',[1/j 1/j  1/j], 'LineWidth',1.7, 'DisplayName',txt2);
 %   %plot3(SX(:,j), SY(:,j), SZ(:,j),'r', 'LineWidth',1.7, 'DisplayName',txt2);
%end

g = surf(SX,SY,SZ,'DisplayName',txt2,'FaceAlpha', 0.8, 'linestyle','none')
colormap('summer')

for i = 1:size(P,1)
    %h(i) = plot3(P(i,:,1),P(i,:,2),P(i,:,3),'r-','LineWidth',1.2,'DisplayName',txt1);
end
for j = 1:size(P,2)
    %plot3(P(:,j,1),P(:,j,2),P(:,j,3),'r-','LineWidth',1.2,'DisplayName',txt1);
end
for j = 1:size(P,2)
    %scatter3(P(:,j,1),P(:,j,2),P(:,j,3),'r','filled');
end
txt0 = ['Derivada'];
bbb = scatter3(DDS(1),DDS(2),DDS(3),'r','filled'); %Derivative
aaa = plot3(derivativeX,derivativeY,derivativeZ,'k-o','LineWidth',1.2,'DisplayName',txt0); %Derivative
%e = 0.05
%for i = 1: length(PX)
%    str = ['P',num2str(i-1)]
%    text(PX(i)+1.5*e,PY(i)+1.5*e,str)
%end
txt3 = ['Pontos de Ruptura'];
%scatter(KX(5:length(U)-k-1), KY(5:length(U)-k-1),'k', 'DisplayName',txt3);
%legend([h(1) g aaa bbb],'Location','north');
%axis([0.8 7.2 1.8 5.2])
hold off