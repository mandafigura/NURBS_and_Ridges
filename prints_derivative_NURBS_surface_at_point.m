clear
format long
%U   = [0 0 0 1 1 1];
%u  = linspace(min(U),max(U),20);
%V   = [0 0 0 0 1 2 3 4 4 4 4];
%v  = linspace(min(V),max(V),150);
%k  = 2;
%l  = 3;

%P          = [0 1 3 3 3 5 5; 0 1 3 3 3 5 5; 0 1 3 3 3 5 5];
%P(:,:,2)   = [0 0 0 0 0 0 0; 3 3 3 3 3 3 3; 6 6 6 6 6 6 6];
%P(:,:,3)   = [2 4 0 0 0 0 4; 2 4 0 0 0 0 4; 2 4 0 0 0 0 4];
%P(:,:,4)   = [1 1 1 1 1 1 1; 1 1 1 1 1 1 1; 1 1 1 1 1 1 1];

%cilinder
%U = [0 0 0 1 1 2 2 3 3 4 4 4]
%u = linspace(min(U),max(U),200);
%V = [0 0 0 1 1 1]
%v  = linspace(min(V),max(V),50);
%k = 2
%l = 2
%p = 1
%q = 1
%point = [0.5 0.5]
%P        = [0 2 4; 0 2 4; 0 2 4; 0 2 4; 0 2 4; 0 2 4; 0 2 4; 0 2 4; 0 2 4]
%P(:,:,2) = [0 0 0; 0 0 0; 2 2 2; 4 4 4; 4 4 4; 4 4 4; 2 2 2; 0 0 0; 0 0 0]
%P(:,:,3) = [2 2 2; 4 4 4; 4 4 4; 4 4 4; 2 2 2; 0 0 0; 0 0 0; 0 0 0; 2 2 2]
%t = sqrt(2)/2;
%P(:,:,4) = [1 1 1; t t t; 1 1 1; t t t; 1 1 1; t t t; 1 1 1; t t t; 1 1 1]

%article cazals
U = [0 0 0 0 0 1 1 1 1 1 ];
u = linspace(min(U),max(U),200);
V = [0 0 0 0 0 1 1 1 1 1 ];
v  = linspace(min(V),max(V),50);
k = 4;
l = 4;
p = 1;
q = 0;
point = [0.8 0.6];

P        = [0 1/4 2/4 3/4 1;   0  1/4 2/4 3/4 4/4;   0  1/4 2/4 3/4 4/4;   0  1/4 2/4 3/4 4/4;  0 1/4 2/4 3/4 1];
P(:,:,2) = [0  0   0   0  0;  1/4 1/4 1/4 1/4 1/4;  2/4 2/4 2/4 2/4 2/4;  3/4 3/4 3/4 3/4 3/4;  1  1   1   1  1];
P(:,:,3) = [0  0   0   0  0;   0   1  -1  -1   0 ;   0  -1   1   1   0 ;   0   1  -1   1   0 ;  0  0   0   0  0];
P(:,:,4) = [1  1   1   1  1;   1   1   1   1   1 ;   1   1   1   1   1 ;   1   1   1   1   1 ;  1  1   1   1  1];


%corrects parameters
P(:,:,1) = transpose(P(:,:,1));
P(:,:,2) = transpose(P(:,:,2));
P(:,:,3) = transpose(P(:,:,3));
P(:,:,4) = transpose(P(:,:,4));


%homogenous coordinates
PW(:,:,1)   = P(:,:,4).*P(:,:,1);
PW(:,:,2)   = P(:,:,4).*P(:,:,2);
PW(:,:,3)   = P(:,:,4).*P(:,:,3);
PW(:,:,4)   = P(:,:,4);

for i = 1:(length(u)-1)
    for j = 1:(length(v)-1)
        [SX(i,j,:)]  = surface_uv_points(k+1,l+1,PW,u(i),v(j),U,V);
    end
end

for i = 1:(length(u)-1)
    for j = 1:(length(v)-1)
        SW(1:4) = SX(i,j,1:4);
        [H] = H_Perspective(SW);
        S(i,j,1:3) = H;
    end
end

%% Derivative
DDS = derivative_order_cd_NURBS_surface(point(1),point(2),U,V,PW,k+1,l+1,1,1);
DD  = derivative_order_cd_NURBS_surface(point(1),point(2),U,V,PW,k+1,l+1,p+1,q+1);

%derivatives (many orders) of the surface
Du   = derivative_order_cd_NURBS_surface(point(1),point(2),U,V,PW,k+1,l+1,2,1);
Dv   = derivative_order_cd_NURBS_surface(point(1),point(2),U,V,PW,k+1,l+1,1,2);
Duu  = derivative_order_cd_NURBS_surface(point(1),point(2),U,V,PW,k+1,l+1,3,1);
Duv  = derivative_order_cd_NURBS_surface(point(1),point(2),U,V,PW,k+1,l+1,2,2);
Dvv  = derivative_order_cd_NURBS_surface(point(1),point(2),U,V,PW,k+1,l+1,1,3);
Duuu = derivative_order_cd_NURBS_surface(point(1),point(2),U,V,PW,k+1,l+1,4,1);
Duuv = derivative_order_cd_NURBS_surface(point(1),point(2),U,V,PW,k+1,l+1,3,2);
Duvv = derivative_order_cd_NURBS_surface(point(1),point(2),U,V,PW,k+1,l+1,2,3);
Dvvv = derivative_order_cd_NURBS_surface(point(1),point(2),U,V,PW,k+1,l+1,1,4);
Duvvv = derivative_order_cd_NURBS_surface(point(1),point(2),U,V,PW,k+1,l+1,2,4);
Dvvvv = derivative_order_cd_NURBS_surface(point(1),point(2),U,V,PW,k+1,l+1,1,5);
Duuvv = derivative_order_cd_NURBS_surface(point(1),point(2),U,V,PW,k+1,l+1,3,3);
Duuuv = derivative_order_cd_NURBS_surface(point(1),point(2),U,V,PW,k+1,l+1,4,2);
Duuuu = derivative_order_cd_NURBS_surface(point(1),point(2),U,V,PW,k+1,l+1,5,1);

%[y]  = sis_func(point(1),point(2));
%[dy] = Dsis_func(point(1),point(2));
%[z]  = sis_func2(point(1),point(2));
%[dz] = Dsis_func2(point(1),point(2));

%first fundamental form
[E,Eu,Ev,F,Fu,Fv,G,Gu,Gv] = first_fundamental_form(Du,Dv,Duu,Duv,Dvv)
%surface normal
[N,barN,barNu,barNv] = surface_normal(Du,Dv,Duu,Duv,Dvv)
%second fundamental form
[e,bare,bareu,barev,f,barf,barfu,barfv,g,barg,bargu,bargv] = second_fundamental_form(Duu,Duv,Dvv,Duuu,Duuv,Duvv,Dvvv,N,barN,barNu,barNv)
%auxiliary terms
[A,Au,Av,B,barB,barBu,barBv,C,barC,barCu,barCv] = auxiliary_ABC(E,Eu,Ev,F,Fu,Fv,G,Gu,Gv,e,bare,bareu,barev,f,barf,barfu,barfv,g,barg,bargu,bargv)
[Pu,Pv,Ru,Rv,Q] = auxiliary_PRQ(A,Au,Av,barB,barBu,barBv,barC,barCu,barCv)
yyy = sis_func(point(1),point(2)) 
NN  = surface_normal(point(1),point(2),U,V,PW,k+1,l+1);

[t1, t2] = principal_directions2d(point(1),point(2),U,V,PW,k+1,l+1)
[T1, T2] = principal_directions3d(point(1),point(2),U,V,PW,k+1,l+1)
[k] = principal_curvatures(point(1),point(2),U,V,PW,k+1,l+1)

%derivada
derivativeX(1) = 0 + DDS(1);
derivativeY(1) = 0 + DDS(2);
derivativeZ(1) = 0 + DDS(3);
exibition   = 1;
derivativeX(2) = exibition.*DD(1) + DDS(1);
derivativeY(2) = exibition.*DD(2) + DDS(2);
derivativeZ(2) = exibition.*DD(3) + DDS(3);

%vetor normal
normalvecX(1) = DDS(1);
normalvecY(1) = DDS(2);
normalvecZ(1) = DDS(3);
normalvecX(2) = N(1) + DDS(1);
normalvecY(2) = N(2) + DDS(2);
normalvecZ(2) = N(3) + DDS(3);

%princ direction blue
princdirec1X(1) = DDS(1);
princdirec1Y(1) = DDS(2);
princdirec1Z(1) = DDS(3);
princdirec1X(2) = T1(1) + DDS(1);
princdirec1Y(2) = T1(2) + DDS(2);
princdirec1Z(2) = T1(3) + DDS(3);

%princ direction red
princdirec2X(1) = DDS(1);
princdirec2Y(1) = DDS(2);
princdirec2Z(1) = DDS(3);
princdirec2X(2) = T2(1) + DDS(1);
princdirec2Y(2) = T2(2) + DDS(2);
princdirec2Z(2) = T2(3) + DDS(3);


hold on
daspect([1 1 1]);
a = [-1.8 -4 1];
[caz,cel] = view(a);
grid on;
%box on;
%daspect([2 5 1]);
xlabel('eixo x');               % legenda no eixo horizontal
ylabel('eixo y');               % legenda no eixo vertical
xlabel('eixo x');               % legenda no eixo horizontal
zlabel('eixo z');               % legenda no eixo vertical
txt1 = ['Pontos de Controle'];
%set(gca,'visible','off')
txt2 = ['Superfície NUBS'];

g = mesh(S(:,:,1),S(:,:,2),S(:,:,3),'DisplayName',txt2);

for i = 1:size(P,1)
    %h(i) = plot3(P(i,:,1),P(i,:,2),P(i,:,3),'k-*','LineWidth',1.2,'DisplayName',txt1);
end
for j = 1:size(P,2)
    %plot3(P(:,j,1),P(:,j,2),P(:,j,3),'k-*','LineWidth',1.2,'DisplayName',txt1);
end
%txt3 = ['Pontos de Ruptura'];
txt0 = ['Derivada'];
bbb = scatter3(DDS(1),DDS(2),DDS(3),'r','filled'); %Derivative
aaa = plot3(derivativeX,derivativeY,derivativeZ,'k-o','LineWidth',1.2,'DisplayName',txt0); %Derivative
%ccc = plot3(normalvecX,normalvecY,normalvecZ,'k-o','LineWidth',1.2,'DisplayName',txt0); %Derivative
ddd = plot3(princdirec1X,princdirec1Y,princdirec1Z,'b-o','LineWidth',1.2,'DisplayName',txt0); %Derivative
eee = plot3(princdirec2X,princdirec2Y,princdirec2Z,'r-o','LineWidth',1.2,'DisplayName',txt0); %Derivativ

%legend([h(1) g aaa],'Location','eastoutside');
hold off
