%ejercicio 1
x = [0 1/3 2/3 1];
Lx = max(x);

fi = inline ("1 + sin(pi/2 * x)");
Nml = inline ("sin(m*pi*x / Lx)","x","Lx","m");
si = inline ("x+1");

opcion = 1; % 1 - Colocacion Puntual / 2 - Galerkin

switch(opcion)
  case 1
    % ------------  Colocacion puntual ----------------------
    %Wl = fi

    %Armamos el sistema Klm*a = f

    Kml = zeros(2,2); %inicializo la matriz
    f = zeros(2,1);%inicializo el vector

    Kml(1,1)=Nml(x(2),Lx,1);
    Kml(1,2)=Nml(x(2),Lx,2);
    Kml(2,1)=Nml(x(3),Lx,1);
    Kml(2,2)=Nml(x(3),Lx,2);

    f(1) = fi(x(2)) - si(x(2));
    f(2) = fi(x(3)) - si(x(3));

    a = inv(Kml) * f;

    fic(1)=fi(x(2));
    fic(2)=fi(x(3));
    fit(1)= si(x(2)) + sum(a(1)*Nml(x(2),Lx,1));
    fit(2)= si(x(3)) + sum(a(2)*Nml(x(3),Lx,2));
  case 2
    
    
    
endswitch

%una ves obtenido el fit se pueden comparar as graficas

figure(1);
plot(fic);
figure(2);
plot(fit);