%ejercicio 2

%-----------DATOS-------------
paso = 1/2;%que tan fina es la subdivision para los caluculos
paso_g = 1/10;%que tan fina es la subdivision para los graficos
b1 = 0;%punto inicial o borde 1
b2 = 1;%punto final o borde 2
%-----------------------------

p = b1:paso:b2;
cant_puntos = length(p);%Cantidad total de puntos segun b1 b2 y paso
mf = 1:cant_puntos;%vector m para la familia de funciones de fourier
X=b1:paso_g:b2;%Vector del dominio mas fino que el calculado para ver que tan buena es la solucion

syms x;

N = x.^(mf+1);
dNdx=diff(N,x);

%Inicializo
Kg=zeros(cant_puntos,cant_puntos);
fg=zeros(cant_puntos,1);

%Galerkin ( Wl = Nl ) debilitando
for l=1:cant_puntos
    for m=1:cant_puntos
         Kg(l,m) = - int(dNdx(l)*dNdx(m),b1,b2) + int(N(l)*N(m),b1,b2) - subs(N(l)*N(m),x,b2) ;
    end
    fg(l) = - int(N(l),b1,b2);
end
ag = Kg\fg;%Solucion de Kg*ag=fg
fi_capa_g = ag'*N'; %Funcion final aproximada con galerkin

%--------Graficos---------
sol_g=subs(fi_capa_g,x,X);%Aproximada con galerkin
plot(sol_g,'r');
legend('Galerkin debilitando');

