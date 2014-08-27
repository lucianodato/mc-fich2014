%ejercicio 1
x = [0 1/3 2/3 1];
Lx = max(x);

%Las funciones inline funcionan poniendo inline ("funcion","argumentos") . Los argumentos se ponen si hay mas de uno sino lo detecta solo
fi = inline ("1 + sin(pi/2 * x)");
Nml = inline ("sin(m*pi*x / Lx)","x","Lx","m");
si = inline ("x+1");

opcion = 2; % 1 - Colocacion Puntual / 2 - Galerkin

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
    %fi original
    fic(1)=fi(x(2));
    fic(2)=fi(x(3));
    %fi capa (aproximado)
    fit(1)= si(x(2)) + sum(a(1:2).*Nml(x(2),Lx,1));
    fit(2)= si(x(3)) + sum(a(1:2).*Nml(x(3),Lx,2));
  case 2
        % ------------  Galerkin ----------------------
    %Wl = Nl

    %Armamos el sistema Klm*a = f
    %Integrales calculadas con el maxima
    
    Kml = zeros(2,2); %inicializo la matriz
    f = zeros(2,1);%inicializo el vector
    
    Kml(1,1)=-(sin(4*pi)-4*pi)/(8*pi);
    Kml(2,2)= -(sin(6*pi)-6*pi)/(12*pi);
    
    f(1) = -(12*pi*sin((5*pi)/2)+15*sin(2*pi)-30*pi*cos(2*pi)-20*pi*sin((3*pi)/2))/(60*pi^2);
    f(2) = -(45*pi*sin((7*pi)/2)+35*sin(3*pi)-105*pi*cos(3*pi)-63*pi*sin((5*pi)/2))/(315*pi^2);
    
    a = inv(Kml)*f;
    
endswitch

%una ves obtenido el fit se pueden comparar as graficas

figure(1);
plot(fic);
figure(2);
plot(fit);