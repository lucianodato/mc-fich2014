function [ temp ] = analitica(cant_celdas,ini,fin,v,Q)
L = fin - ini;
h = L/cant_celdas;%distancia que separa las celdas (si son regulares)

%Solver analitico. By RJM. 
syms y(t)
Dy = diff(y);
Dy_2 = diff(y,2);
%Aca remplazmos con la ecuacion.
s = dsolve(v*Dy - Dy_2  - Q  == 0, y(0) == 0, Dy(1) == 1);
temp = subs(s,t,0:h:L-h);
end

