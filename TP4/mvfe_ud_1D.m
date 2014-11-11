%Resolucion estacionaria por metodo de volumen finito para la ecuacion de
%transmision de Calor en 1D con termino advectivo por upwind difference y
%con celdas de tamaño regular

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
cm_h = 1;%h de la condicion mixta
cm_k = 1;%k de la condicion mixta si la hay
cm_finf = 1;%temperatura externa fi inf

%Definicion de las condiciones de borde (-1 significa que no aplica)
cbd_i = 0;%condicion de borde dirichlet izquierda
cbd_d = 1;%condicion de borde dirichlet derecha
cbn_i = -1;%condicion de borde neumann izquierda
cbn_d = -1;%condicion de borde neumann derecha
cbm_i = -1;%condicion de borde mixta izquierda
cbm_d = -1;%condicion de borde mixta derecha

%Matriz del lado derecho (todos los terminos que tienen incognitas fi)
A = zeros(cant_celdas,cant_celdas);

for i = 1:cant_celdas
    %Partes de la ecuacion de transmision que colaboran con la matriz
    switch i
        case 1
            %Caso que es la primera celda
            if (cbd_i ~= -1) %Condicion Dirichlet
                A(i,i) = -v + k/h * -3;
                A(i,i+1) = k/h;
            else
                if(cbn_i ~= -1)
                    %Condicion Neumann
                    A(i,i) = -2*v + k/h * -3;
                    A(i,i+1) = k/h;
                else
                    %Condicion Mixta
                    A(i,i) = -v * 1/2 + k/h * -3;
                    A(i,i+1) = -v * 1/2 + k/h;
                end
            end
        case cant_celdas
            %Caso ultima celda
            if (cbd_d ~= -1) %Condicion Dirichlet
                A(i,i-1) = -v + k/h;
                A(i,i) = -v + k/h * -3;
            else
                if(cbn_d ~= -1)
                    %Condicion Neumann
                    A(i,i-1) = -v + k/h;
                    A(i,i) = -v + k/h * -3;
                else
                    %Condicion Mixta
                    A(i,i) = -v * 1/2 + k/h * -3;
                    A(i,i+1) = -v * 1/2 + k/h;
                end
            end
        otherwise
            %celdas interna
            A(i,i-1) = - v + k/h;
            A(i,i) = -v + k/h * -2;
            A(i,i+1) = k/h;
    end
end


%vector del lado izquierdo (todos los terminos que no tienen icognita fi y condiciones de borde)
b = zeros(cant_celdas,1);

for i = 1:cant_celdas
    %Partes de la ecuacion de transmision que colaboran con la matriz
    switch i
        case 1
            %Caso que es la primera celda
            if (cbd_i ~= -1) 
                %Condicion Dirichlet
                b(i) = -Q*h - 2*k/h * cbd_i + v * cbd_i;
            else
                if (cbn_i ~= -1)
                    %Condicion Neumann
                    b(i) = -Q*h + k*cbn_i + v * (-1*cbn_i);%el termino advectivo actua en la direccion de la cara
                else
                    %Condicion Mixta
                    b(i) = -Q*h - 2*k/h * cbn_i + cm_h*cm_finf;
                end
            end
        case cant_celdas
            %Caso ultima celda
            if (cbd_d ~= -1)
                %Condicion Dirichlet
                b(i) = -Q*h - 2*k/h * cbd_d + v * cbd_d;
            else
                if (cbn_d ~= -1)
                    %Condicion Neumann
                    b(i) = -Q*h + k*cbn_d ;%+ v * (1*cbn_d);%el termino advectivo actua en la direccion de la cara
                else
                    %Condicion Mixta
                    b(i) = -Q*h - 2*k/h * cbn_d ;
                end
            end
        otherwise
            %celdas interna
            b(i) = -Q*h;
    end
end

%Resolucion del sistema
temp=A\b;
