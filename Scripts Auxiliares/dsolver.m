%Solver analitico. By RJM. 
syms y(t)
Dy = diff(y);
Dy_2 = diff(y,2);
%Acá remplazmos con la ecuación.
s = dsolve(100*Dy - Dy_2  - 1  == 0, y(0) == 0, Dy(1) == 1);
%Ploteamos la func.
hold on;
plot(subs(s,t,0:0.2:1),'g');

