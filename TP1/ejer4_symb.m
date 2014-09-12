%ejercicio 3

%-----------DATOS-------------
paso = 1/2;%que tan fina es la subdivision para los caluculos
paso_g = 1/10;%que tan fina es la subdivision para los graficos
nu=1/4;
E=30/16;
a1 = 1;%Lado X derecha
a2 = -1;%Lado X izquierda
b1 = 1;%Lado Y derecha
b2 = -1;%Lado Y izquierda
%-----------------------------

cant_puntos = 3;%Cantidad total de puntos 
mf = 1:cant_puntos;%vector m para la familia de funciones de fourier
p = a1:paso_g:a2;%Eje X
q = b1:paso_g:b2;%Eje Y
[X,Y]=meshgrid(p,q);%Vectores del dominio mas fino que el calculado para ver que tan buena es la solucion

syms u v x y;

tx=E*(1-y^2)/(1+nu);
ty=0;
D=(E/(1-nu^2))*[1 nu 0; nu 1 0; 0 0 (1-nu)/2]; %Tension plana
t=[tx ty]';
N11=x*(1-y^2);
N21=(x^2)*(1-y^2);
N31=x^3*(1-y^2);
N12=y*(1-y^2);
N22=(y^2)*(1-y^2);
N32=(y^3)*(1-y^2);
N=[N11 N12; N21 N22; N31 N32];
N2=[N11 N21 N31 N12 N22 N32];

%Inicializo
Kg=zeros(cant_puntos,cant_puntos);
fg=zeros(2*cant_puntos,1);

%Galerkin ( Wl = Nl )
%Esto es la ecuacion 1 y 2 y al final se arma la matriz K que es simetrica por eso solo se calcula la ecuacion 1
for l=1:cant_puntos
    for m=1:cant_puntos
        LNl=[diff(N(l,1),x,1) 0; 0 diff(N(l,2),y,1); diff(N(l,1),y,1) diff(N(l,2),x,1)];
        LNm=[diff(N(m,1),x,1) 0; 0 diff(N(m,2),y,1); diff(N(m,1),y,1) diff(N(m,2),x,1)];
        Klm=int(int(LNl'*D*LNm,x,a2,a1),y,b2,b1);
        Kg(2*l-1:2*l,2*m-1:2*m)=Klm;
    end
    %Para la primera ecuacion vale para la segunda es 0
    if (l <= cant_puntos)
        fg(l) = 6 * int((1-y^2)*subs(N(l,1),x,1),y,b2,b1);
    end
end
ag = Kg\fg;%Solucion de Kg*ag=fg
