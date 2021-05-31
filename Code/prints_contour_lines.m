load('100seed100')
tamanho = 100;
px = linspace (0,1, tamanho)';
py = linspace (0,1, tamanho)';
[x, y] = meshgrid (px, py);
z1 = zeros(tamanho);
z2 = zeros(tamanho);
for i=1:tamanho-1
    for j=1:tamanho-1
        [z1(i,j) z2(i,j)] = principal_curvatures([x(i,j) y(i,j)]);
    end
end

figure(1)
hold on
daspect([1 1 1]);
contour(x,y,z1,50,'blue')
scatter(testk1(:,1),testk1(:,2),'k')
axis([0 1 0 1])
hold off

figure(2)
hold on
daspect([1 1 1]);
contour(x,y,z2,50,'red')
scatter(testk2(:,1),testk2(:,2),'k')
axis([0 1 0 1])
hold off