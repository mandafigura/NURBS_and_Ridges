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

%U = [0 0 0 1 1 2 2 3 3 4 4 4]
%u = linspace(min(U),max(U),200);
%
%V = [0 0 0 1 1 1]
%v  = linspace(min(V),max(V),50);
%k = 2
%l = 2

%P        = [0 2 4; 0 2 4; 0 2 4; 0 2 4; 0 2 4; 0 2 4; 0 2 4; 0 2 4; 0 2 4]
%P(:,:,2) = [0 0 0; 0 0 0; 2 2 2; 4 4 4; 4 4 4; 4 4 4; 2 2 2; 0 0 0; 0 0 0]
%P(:,:,3) = [2 2 2; 4 4 4; 4 4 4; 4 4 4; 2 2 2; 0 0 0; 0 0 0; 0 0 0; 2 2 2]
%t = sqrt(2)/2;
%P(:,:,4) = [1 1 1; t t t; 1 1 1; t t t; 1 1 1; t t t; 1 1 1; t t t; 1 1 1]


%PW(:,:,1)   = P(:,:,4).*P(:,:,1);
%PW(:,:,2)   = P(:,:,4).*P(:,:,2);
%PW(:,:,3)   = P(:,:,4).*P(:,:,3);
%PW(:,:,4)   = P(:,:,4);

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
        [SX(i,j,:)]  = Ponto_uv_da_superficie_nurbs(k+1,l+1,PW,u(i),v(j),U,V);
    end
end

for i = 1:(length(u)-1)
    for j = 1:(length(v)-1)
        SW(1:4) = SX(i,j,1:4);
        [H] = H_Perspectiva(SW);
        S(i,j,1:3) = H;
    end
end

DDS = derivative_order_cd_NURBS_surface(0.1,0.2,U,V,PW,k+1,l+1,1,1);
[T1,T2] = principal_directions3d(0.1,0.2)

T1 = T1/norm(T1)
T2 = T2/norm(T2)

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
%daspect([1 1 1]);
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

g = mesh(S(:,:,1),S(:,:,2),S(:,:,3),'DisplayName',txt2)
ddd = plot3(princdirec1X,princdirec1Y,princdirec1Z,'b-o','LineWidth',1.2); %Derivative
eee = plot3(princdirec2X,princdirec2Y,princdirec2Z,'r-o','LineWidth',1.2); %Derivativ

txt3 = ['Pontos de Ruptura'];
legend([g],'Location','eastoutside');
hold off

