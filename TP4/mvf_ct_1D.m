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

%Matriz del lado derecho (todos los terminos que tienen incognitas fi)
A = zeros(cant_celdas,cant_celdas);

for i = 1:cant_celdas
    %Partes de la ecuacion de transmision que colaboran con la matriz
    switch i
        case 1
            %Caso que es la primera celda
            if (cbd_i ~= -1) %Condicion Dirichlet
                A(i,i) = v * 1/2 + k/h * -3;
                A(i,i+1) = v * 1/2 + k/h;
            else %Condicion Neumann
                A(i,i) = v + v * 1/2 + k/h * -3;
                A(i,i+1) = v * 1/2 + k/h;
            end
        case cant_celdas
            %Caso ultima celda
            if (cbd_d ~= -1) %Condicion Dirichlet
                A(i,i-1) = v * 1/2 + k/h;
                A(i,i) = v * 1/2 + k/h * -3;
            else %Condicion Neumann
                A(i,i-1) = v * 1/2 + k/h;
                A(i,i) = v + v * 1/2 + k/h * -3;
            end
        otherwise
            %celdas interna
            A(i,i-1) = v * 1/2 + k/h;
            A(i,i) = v + k/h * -2;
            A(i,i+1) = v * 1/2 + k/h;
    end
end


%vector del lado izquierdo (todos los terminos que no tienen icognita fi y condiciones de borde)
b = zeros(cant_celdas,1);

for i = 1:cant_celdas
    %Partes de la ecuacion de transmision que colaboran con la matriz
    switch i
        case 1
            %Caso que es la primera celda
            if (cbd_i ~= -1) %Condicion Dirichlet
                b(i) = -Q*h - 2*k/h * cbd_i - v * cbd_i;
            else %Condicion Neumann
                b(i) = -Q*h - 2*k/h * cbn_i - v * (-1*cbn_i);%el termino advectivo actua en la direccion de la cara
            end
        case cant_celdas
            %Caso ultima celda
            if (cbd_d ~= -1)
                b(i) = -Q*h - 2*k/h * cbd_d;
            else
                b(i) = -Q*h - 2*k/h * cbn_d - v * (1*cbn_i);%el termino advectivo actua en la direccion de la cara
            end
        otherwise
            %celdas interna
            b(i) = -Q*h;
    end
end

%Resolucion del sistema
temp=A\b;
