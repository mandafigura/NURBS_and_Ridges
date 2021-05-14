load('100seed100')
c1 = principal_curvature_k1(testk2(2,1),testk2(2,2))
tamanho = 100
epsilon = 0.2
%px = linspace ((c1-epsilon), (c1+epsilon), tamanho)';
%py = linspace ((c1-epsilon), (c1+epsilon), tamanho)';
px = linspace (0,1, tamanho)';
py = linspace (0,1, tamanho)';
[x, y] = meshgrid (px, py);
z = zeros(tamanho);
for i=1:tamanho-1
    for j=1:tamanho-1
        z(i,j) = principal_curvature_k1(x(i,j),y(i,j));
    end
end

hold on
scatter(testk2(:,1),testk2(:,2),'r')
contour(x,y,z,50)
axis([0 1 0 1])
hold off