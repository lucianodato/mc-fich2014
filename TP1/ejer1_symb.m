%ejercicio 1

%-----------DATOS-------------
paso = 1/3;%que tan fina es la subdivision para los caluculos
paso_g = 1/10;%que tan fina es la subdivision para los graficos
b1 = 0;%punto inicial o borde 1
b2 = 1;%punto final o borde 2
%-----------------------------

p = b1:paso:b2;
L = max(p);%Maximo punto para la normalizacion de furier
cant_puntos = length(p);%Cantidad total de puntos segun b1 b2 y paso
mf = 1:cant_puntos;%vector m para la familia de funciones de fourier
X=b1:paso_g:b2;%Vector del dominio mas fino que el calculado para ver que tan buena es la solucion

syms x;

fi = 1 + sin((pi/2) * x);
N = sin(mf.*pi*x/ L);%Familia de Fourier
psi = (x+1);

%Inicializo
Kg=zeros(cant_puntos,cant_puntos);
fg=zeros(cant_puntos,1);
Kcp=zeros(cant_puntos,cant_puntos);
fcp=zeros(cant_puntos,1);

%Buscamos los pesos a para cada caso

%Galerkin ( Wl = Nl )
for l=1:cant_puntos
    for m=1:cant_puntos
         Kg(l,m) = int(N(l)*N(m),b1,b2); 
    end
    fg(l) = int(N(l)*(fi-psi),b1,b2);
end
ag = Kg\fg;%Solucion de Kg*ag=fg
fi_capa_g = psi + ag'*N'; %Funcion final aproximada con galerkin

%Colocacion Puntual ( Wl = f(x(l)) )
for l=1:cant_puntos
    for m=1:cant_puntos
         Kcp(l,m) = subs(N(m),x,p(l)); 
    end
    fcp(l) = subs((fi-psi),x,p(l)); 
end
acp = Kcp\fcp;%Solucion de Kcp*acp=fcp
fi_capa_cp = psi + acp'*N';%Funcion final aproximada con colocacion puntual


%--------Graficos---------
sol_g=subs(fi_capa_g,x,X);%Aproximada con galerkin
sol_cp=subs(fi_capa_cp,x,X);%Aproximada con colocacion puntual
fun=subs(fi,x,X);%Original
hold on;
plot(sol_g,'r');
plot(sol_cp,'g');
plot(fun,'b');
legend('Galerkin','Colocacion Puntual','Funcion Original');
