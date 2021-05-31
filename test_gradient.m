load('100seed100')
[m1,n1] = size(testk1);
[m2,n2] = size(testk2);
bla = [0 0]
hold on
for i=1:m1
   [gradient] = ridge_gradient(testk1(i,:),'blue');
   gradient = gradient./norm(gradient);
   bla(1,:) = testk1(i,:);
   bla(2,:) = ((0.1).*gradient + testk1(i,:));
   scatter(bla(1,1),bla(1,2),'b')
   plot(bla(:,1),bla(:,2),'b')
end
daspect([1 1 1]);
axis([0 1 0 1])
hold off

hold on
for i=1:m2
   [gradient] = ridge_gradient(testk2(i,:),'red');
   gradient = gradient./norm(gradient);
   bla(1,:) = testk2(i,:);
   bla(2,:) = ((0.1).*gradient + testk2(i,:));
   scatter(bla(1,1),bla(1,2),'r')
   plot(bla(:,1),bla(:,2),'r')
end
daspect([1 1 1]);
axis([0 1 0 1])
hold off

