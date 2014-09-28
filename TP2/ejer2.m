syms x Vh a b c x1 x2 x3 u1 u2 u3 h;
%funcion generica
Vh=a*x^2+b*x+c;
%divido el elemento en 3 puntos para armar el sistema de ecuaciones
x1=-h/2;
x2=0*h;
x3=h/2;
%armo la matriz considerando cada renglon como la ecuacion de cada punto
M = [x1^2 x1 1;x2^2 x2 1;x3^2 x3 1];
%armo el vector solucion considerando lo correcto para cada punto
%que quede el punto en 1 y el resto en 0
u1 = [1 0 0]';
u2 =[0 1 0]';
u3 = [0 0 1]';
e1 = M\u1;
e2 = M\u2;
e3 = M\u3;

N1 = subs(Vh,[a,b,c],[e1(1),e1(2),e1(3)]);
N2 = subs(Vh,[a,b,c],[e2(1),e2(2),e2(3)]);
N3 = subs(Vh,[a,b,c],[e3(1),e3(2),e3(3)]);