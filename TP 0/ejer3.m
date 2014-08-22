function [S,E,x] = estadoplano (Lx,Ly,C)

	//Declaración de variables,operador diferencial L y syms.
	//Lx y Ly las recibo como parametros.
	// S es el tensor de tensiones y E el tensor de deformación.
	//El valor de C es la constante de Young del material. La recibo como parametro.

	u = syms f(x, y);
	v= syms f(x, y);

	u(x, y)=-0.1*sin(pi/Lx).cos(pi*x/2*Ly);
	v(x, y)=-0.1*sin(pi/2*Lx).cos(pi*x/Ly);

	x= zeros (1,2); // fuerzas de cuerpo Xi

	E= zeros(2,2);
	S=E;

	//------__-------

	//(A)- CALCULO EL TENSOR DE DEFORMACIÓN

	Ex=diff(u,'x');
	Ey=diff(v,'y');
	Exy=diff(u,'y') + diff(v,'x');
	Eyx= Exy;

	//GRAFICO EL TENSOR DE DEFORMACIÓN. NO SE COMO SE HACE... JE
		
	figure(1);

		//IMAGINO Q TENGO Q ASIGNARLE VALORES A 'X' E 'Y',ARMAR LA MATRIZ Y GRAFICAR DE ALGUNA MANERA.
	plot();

	//------__-------
	
	//(B)- CALCULO EL TENSOR DE TENSIONES.

	Sx=C*.Ex;
	Sy=C*Ey;
	Sxy= C*Exy;
	Syx= Sxy;

	//GRAFICO EL TENSOR DE TENSIONES.
	
	figure(2);

		//IMAGINO Q TENGO Q ASIGNARLE VALORES A 'X' E 'Y',ARMAR LA MATRIZ Y GRAFICAR DE ALGUNA MANERA.
	plot();

	//------__-------

	//(C)- CALCULO LAS FUERZAS DE CUERPO.

	dSx=diff(Sx,'x');
	dSy=diff(Sy,'y');
	dSxy= diff(Sxy,'y') + diff(Syx,'x');
	
	X(1) = - dSx - dSxy ; //Resp. a x
	X(2) = - dSy - dSyx; // Resp. a y

ENDFUNCTION

