%Solver analitico. By RJM. 
syms y(t)
Dy = diff(y);
%Ac� remplazmos con la ecuaci�n.
s = dsolve(diff(y, 2)  - 10 == 0, y(0) == 1, Dy(1) == 1);
%Ploteamos la func.
plot(subs(s,t,0:0.1:1));