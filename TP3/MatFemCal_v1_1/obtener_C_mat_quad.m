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
        syms chi nu;

=======
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
>>>>>>> de61ad3c078db6a5c085eee4e37eff9d5e4a3c2a
        xi = nodes(1,1);
        yi = nodes(1,2);
        xj = nodes(2,1);
        yj = nodes(2,2);
        xk = nodes(3,1);
        yk = nodes(3,2);
        xl = nodes(4,1);
        yl = nodes(4,2);
>>>>>>> origin/master
        
        %Familia de forma cuadrangular
        Ni = 1/4 * (1-chi)*(1-nu);%symbolic
        Nj = 1/4 * (1+chi)*(1-nu);
        Nk = 1/4 * (1-chi)*(1-nu);
        Nl = 1/4 * (1-chi)*(1+nu);    
        
        N = [Ni Nj Nk Nl];
        xx = xi*Ni + xj*Nj + xk*Nk + xl*Nl; 
        yy = yi*Ni + yj*Nj + yk*Nk + yl*Nl;
        
        dxxchi = diff(xx,chi);
        dxxnu = diff(xx,nu);
        dyychi = diff(yy,chi);
        dyynu = diff(yy,nu);
        
        J = [dxxchi dyychi;dxxnu dyynu];
        detJ = det(J);
        
        M = [];
        for i=1:4
            for j=1:4
                M(i,j) = int(int(N(i)*N(j)*detJ,nu,-1,1),chi,-1,1);
            end
        end
        
        Elem_C_Mat = cp*rho*M;
end