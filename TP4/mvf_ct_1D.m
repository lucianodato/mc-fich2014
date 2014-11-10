%Definicion de parametros globales
cant_celdas = 5; %numeros de elementos en la recta 1D
ini = 0;%punto de inicio de la recta
fin = 1;%punto de fin de la recta
L = fin - ini;
h = L/cant_celdas;%distancia que separa las celdas (si son regulares)

%Definicion de parametros especificos de la ecuacion de balance termico
k = 1;%Constante de difusividad
Q = 1;%fuente
v = 1;%velocidad
c = 1;%constante reactiva

%Definicion de las condiciones de borde (-1 significa que no aplica)
cbd_i = 1;%condicion de borde dirichlet izquierda
cbd_d = -1;%condicion de borde dirichlet derecha
cbn_i = -1;%condicion de borde neumann izquierda
cbn_d = 1;%condicion de borde neumann derecha

%Problema estacionario
A = zeros(cant_celdas,cant_celdas);
b = zeros(cant_celdas,1);

%Matriz del lado derecho (todos los terminos que tienen incognitas fi)
for i = 1:cant_celdas
    %Partes de la ecuacion de transmision que colaboran con la matriz
    switch cant_celdas
        case 1
            %Caso que es la primera celda
            A(i,i) = k/h * -3;
            A(i,i+1) = k/h;
        case cant_celdas
            %Caso ultima celda
            A(i,i-1) = k/h;
            A(i,i) = k/h * -3;
        otherwise
            %celdas interna
            A(i,i-1) = k/h;
            A(i,i) = k/h * -2;
            A(i,i+1) = k/h;
    end
end
%vector del lado izquierdo (todos los terminos que no tienen icognita fi y condiciones de borde)

for i = 1:cant_celdas
    %Partes de la ecuacion de transmision que colaboran con la matriz
    switch cant_celdas
        case 1
            %Caso que es la primera celda
            if (cbd_i ~= -1)
                b(i) = -Q*h - 2*k/h * cbd_i;
            else
                b(i) = -Q*h - 2*k/h * cbn_i;
            end
        case cant_celdas
            %Caso ultima celda
                        if (cbd_i ~= -1)
                b(i) = -Q*h - 2*k/h * cbd_d;
            else
                b(i) = -Q*h - 2*k/h * cbn_d;
            end
        otherwise
            %celdas interna
            b(i) = -Q*h;
    end
end
