%ejercicio 3

%-----------DATOS-------------
paso = 1/2;%que tan fina es la subdivision para los caluculos
paso_g = 1/10;%que tan fina es la subdivision para los graficos
a1 = 1;%Lado X derecha
a2 = -1;%Lado X izquierda
b1 = 1;%Lado Y derecha
b2 = -1;%Lado Y izquierda
%-----------------------------

cant_puntos = abs((a2-a1)/paso) +1;%Cantidad total de puntos segun b1 b2 y paso
mf = 1:cant_puntos;%vector m para la familia de funciones de fourier
p = a1:paso_g:a2;%Eje X
q = b1:paso_g:b2;%Eje Y
[X,Y]=meshgrid(p,q);%Vectores del dominio mas fino que el calculado para ver que tan buena es la solucion

syms x y;

psi = (1-x^2) + (1-y^2);%psi que cumple con las condiciones de borde necesarias
N=sin(pi.*mf*x/(a1-a2)) .* sin(pi.*mf*y/(b1-b2));%familia de funciones que cumple con las condiciones de borde
d2Ndx2=diff(N,x,2);
d2Ndy2=diff(N,y,2);

%Inicializo
Kg=zeros(cant_puntos,cant_puntos);
fg=zeros(cant_puntos,1);

%Galerkin ( Wl = Nl )
for l=1:cant_puntos
    for m=1:cant_puntos
         Kg(l,m) = int(int(N(l)*(d2Ndx2(m)+d2Ndy2(m)),x,a2,a1),y,b2,b1);
    end
    fg(l) = 4 * int(int(N(l),x,a2,a1),y,b2,b1);
end
ag = Kg\fg;%Solucion de Kg*ag=fg
fi_capa_g = psi + ag'*N'; %Funcion final aproximada con galerkin

%--------Graficos---------
sol_g=subs(fi_capa_g,[x,y],[X,Y]);
figure(1);
mesh(X,Y,sol_g);
legend('Galerkin con a calculados');
figure(2);
mesh(X,Y,subs(psi,[x,y],[X,Y]));
legend('Galerkin con a todos ceros');