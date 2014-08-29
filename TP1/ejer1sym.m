%ejercicio 1
p = [0 1/3 2/3 1];
Lp = max(p);

syms m l x;

fi = 1 + sin(pi/2 * x);
Nm = sin(m*pi*x / Lp);%Familia de Fourier
Nl = sin(l*pi*x / Lp);
si = x+1;

Kml = zeros(2,2); %inicializo la matriz
f = zeros(1,2);%inicializo el vector

%Buscamos los pesos a para cada caso

opcion = 2; % 1 - Colocacion Puntual / 2 - Galerkin

switch(opcion)
  case 1
    % ------------  Colocacion puntual ----------------------
    %Wl = fi
    %En este caso no necesitamos integrar solo evaluar la funcion en el punto
    %Consideramos solo el dominio interior de omega ya que las condiciones de borde estan cubiertas por psi
    %Armamos el sistema Klm*a = f

    Kml(1,1)=subs(Nm,[m,x],[1,p(2)]);
    Kml(1,2)=subs(Nm,[m,x],[2,p(2)]);
    Kml(2,1)=subs(Nm,[m,x],[1,p(3)]);
    Kml(2,2)=subs(Nm,[m,x],[2,p(3)]);

    f(1) = subs(fi,x,p(2)) - subs(si,x,p(2));
    f(2) = subs(fi,x,p(3)) - subs(si,x,p(3));
    
    a = f/Kml;
  case 2
    % ------------  Galerkin ----------------------
    %Wl = Nl

    %Armamos el sistema Klm*a = f
    %Primero integramos y hallamos la forma de la matriz Klm y del vector respuesta f
    
    %Aunque para a familia de funciones de furier el metodo de galerkin ya esta definido en una ecuacion, calculamos de todas maneras las integrales
    
    Nlm = Nl * Nm;
    Nlmi = int(Nlm);%integral indefinida porque despues reemplazo los valores en el calculo de cada celda de la matriz 
    
    fl = Nl * (fi - si);
    fli = int(fl);
    
    Kml(1,1)=subs(Nlmi,[m,l,x],[1,1,p(2)]);
    Kml(1,2)=subs(Nlmi,[m,l,x],[1,2,p(2)]);
    Kml(2,2)=subs(Nlmi,[m,l,x],[2,2,p(3)]);
    Kml(2,1)=subs(Nlmi,[m,l,x],[2,1,p(3)]);
    
    f(1) = subs(fli,[l,x],[1,p(2)]);
    f(2) = subs(fli,[l,x],[2,p(3)]);
    
    a = f/Kml;
    
end

%una ves obtenido el fit se pueden comparar as graficas
m=[1,2,3,4];
an = [1,a(1),a(2),2]
fic = 1 + sin(pi/2 .* p);
fit = p + 1 + sum(an.*sin(pi/Lp * m.*p));

figure(1);
plot(fic);
hold on;
plot(fit);