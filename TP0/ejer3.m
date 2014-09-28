
%Declaración de variables,operador diferencial L y syms.
%Lx y Ly las recibo como parametros.
%S es el tensor de tensiones y E el tensor de deformación.
%El valor de C es la constante de Young del material. La recibo como parametro.


%Bordes placa placa
Lx=1;
Ly=1;

%Declaro las variables simbolicas
symbols
x=sym("x");
y=sym("y");

%Desplazamientos U(u,v)
u=-0.1*Sin(pi/Lx)*Cos(pi*x/2*Ly);
v=-0.1*Sin(pi/2*Lx)*Cos(pi*x/Ly);

po = 0.3; % poisson coef. para el metal
yo=210000; %modulo de young

coe = yo / (1-po^2); % factor q multiplica  a los miembros de la matriz D

D = zeros(3,3);
D(1,1) =1;
D(1,2) = po;
D(2,1) = po;
D(2,2) = 1;
D(3,3) = (1-po)/2;

D= (yo/(1-po^2))*D; % matriz D relación tensión - deformación para elasticidad

%------__-------

%(A)- CALCULO EL TENSOR DE DEFORMACIÓN

Ex=differentiate(u,x);
Ey=differentiate(v,y);
Eyx=differentiate(u,y);
Exy=Eyx;


%GRAFICO EL TENSOR DE DEFORMACIÓN. NO SE COMO SE HACE... JE
	
%figure(1);

	%IMAGINO Q TENGO Q ASIGNARLE VALORES A 'X' E 'Y',ARMAR LA MATRIZ Y GRAFICAR DE ALGUNA MANERA.
%plot();

%------__-------

%(B)- CALCULO EL TENSOR DE TENSIONES.
Sx=coe*Ex;
Sy=coe*(po*Ex + po*Ey);

aux = (1-po)/2;
Syx= coe*(aux*Eyx + aux*Exy);
Sxy = Syx;

%GRAFICO EL TENSOR DE TENSIONES.

%figure(2);

	%IMAGINO Q TENGO Q ASIGNARLE VALORES A 'X' E 'Y',ARMAR LA MATRIZ Y GRAFICAR DE ALGUNA MANERA.
%plot();

%------__-------

%(C)- CALCULO LAS FUERZAS DE CUERPO.

dSx=differentiate(Sx,x);
dSy=differentiate(Sy,y);
dSyx= differentiate(Syx,y) + differentiate(Sxy,x);
dSxy=dSyx;

bx = -1*(dSx + dSxy) ; %Resp. a x
by = -1*(dSy + dSyx); % Resp. a y

