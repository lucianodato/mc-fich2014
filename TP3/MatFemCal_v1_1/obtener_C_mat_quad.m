% Las entradas vienen en este formato
%
% coord =
% 
%    40.0000   54.2857
%    40.0000   51.4286
%    42.5981   52.8571
% 

function Elem_C_Mat = obtener_C_mat_quad(nodes, rho, cp)

<<<<<<< HEAD
fform = @(s,t)[(1-s-t+s*t)/4,(1+s-t-s*t)/4,(1+s+t+s*t)/4,(1-s+t-s*t)/4];

pospg = [ -0.577350269189626E+00 , 0.577350269189626E+00 ];
pespg = [  1.0E+00 , 1.0E+00]; % Esta variable es el peso de wi y wj, pesos relativos a nu y chi en la integración
M = zeros(4,4);

for i=1 : 2
    for j=1 : 2
        lcffm = fform(pospg(i),pospg(j)) ;    % FF at gauss point
        xjacm = lcffm*nodes ;                 % Jacobian matrix
        ctder = xjacm\lcffm ;                 % FF Cartesian derivates
        darea = det(xjacm)*pespg(i)*pespg(j);
=======
        %Area del cuadrangulo (con las cordenadas de cada nodo de elemento)
        xi = nodes(1,1);
        yi = nodes(1,2);
        xj = nodes(2,1);
        yj = nodes(2,2);
        xk = nodes(3,1);
        yk = nodes(3,2);
        xl = nodes(4,1);
        yl = nodes(4,2);
>>>>>>> origin/master
        
        l1= 1/2*sqrt((xj-xi)^2 + (yj-yi)^2);
        l2= 1/2*sqrt((xk-xj)^2 + (yk-yj)^2);
        A = l1*l2;
        
        %Familia de forma cuadrangular
        Ni = ((l1-x)*(l2-y))/(4*l1*l2);%symbolic
        Nj = ((l1+x)*(l2-y))/(4*l1*l2);
        Nk = ((l1+x)*(l2+y))/(4*l1*l2);
        Nl = ((l1-x)*(l2+y))/(4*l1*l2);
        
        %For every side there is one matrix

        


        Elem_C_Mat = rho * cp * M;

end