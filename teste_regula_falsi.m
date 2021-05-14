clear
format long
tic
[x,y,z,flag,u] = check_simplexes([0.42 0.62],[0.47 0.62],[0.42 0.6],[0.47 0.6])
[answer] = regula_falsi_2d(x,y,z,3);
[bbb n status] = Newton_Sis([0.47 0.60],1.0e-5,500);

toc
if ((0.333 <= 1.001) && (4 <= 1.0001))
    fprintf('\n deucerto\n')
else
    fprintf('\n NÃO deucerto\n')
end