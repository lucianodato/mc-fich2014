%Solver analitico. By RJM. 
syms y(t)
Dy = diff(y);
%Ac� remplazmos con la ecuaci�n.
s = dsolve(100*diff(y, 1) - diff(y, 2)  - 1 == 0, y(0) == 0, Dy(1) == 1);
%Ploteamos la func.
hold on;
plot(subs(s,t,0:0.2:0.8),'g');