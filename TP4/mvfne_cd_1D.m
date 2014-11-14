%Resolucion no estacionaria por metodo de volumen finito para la ecuacion de
%transmision de Calor en 1D con termino advectivo por central difference y
%con celdas de tamaño regular y utilizando crank-nicholson

%Definicion de parametros globales
cant_celdas = 5; %numeros de elementos en la recta 1D
ini = 0;%punto de inicio de la recta
fin = 1;%punto de fin de la recta
L = fin - ini;
h = L/cant_celdas;%distancia que separa las celdas (si son regulares)

%Definicion de parametros especificos de la ecuacion de balance termico
k = -1;%Constante de difusividad
Q = -1;%fuente
v = -100;%velocidad
cm_h = 1;%h de la condicion mixta
cm_k = 1;%k de la condicion mixta si la hay
cm_finf = 1;%temperatura externa fi inf

%Definicion de parametros no estacionarios
dt = 2*10^-3;%paso de tiempo
t_max = 0.01;%tiempo maximo
t_ini = 0;%tiempo inicial
cant_pasos_tiempo = (t_max-t_ini)/dt;
temp_t = zeros(cant_celdas,1); %Por el inicial

%Definicion de las condiciones de borde (-1 significa que no aplica)
cbd_i = 0;%condicion de borde dirichlet izquierda
cbd_d = -1;%condicion de borde dirichlet derecha
cbn_i = -1;%condicion de borde neumann izquierda
cbn_d = 1;%condicion de borde neumann derecha
cbm_i = -1;%condicion de borde mixta izquierda %cualquier numero distinto de -1 la activa
cbm_d = -1;%condicion de borde mixta derecha



for i = 1:cant_pasos_tiempo
    
    %Matriz del lado derecho (todos los terminos que tienen incognitas fi)
    A = zeros(cant_celdas,cant_celdas);
    
    for j = 1:cant_celdas
        %Partes de la ecuacion de transmision que colaboran con la matriz
        switch j
            case 1
                %Caso que es la primera celda
                if (cbd_i ~= -1) %Condicion Dirichlet
                    A(j,j) = h/dt + 1/2*(-v * 1/2 + k/h * -3);
                    A(j,j+1) = 1/2 * (-v * 1/2 + k/h);
                else
                    if(cbn_i ~= -1)
                        %Condicion Neumann
                        A(j,j) = h/dt +1/2*(v - v * 1/2 + k/h);
                        A(j,j+1) = 1/2*(-v * 1/2 + k/h);
                    else
                        %Condicion Mixta
                        difer=4*cm_k/(2*h*cm_k+cm_h*h^2) - 2/h;
                        A(j,j) = h/dt + 1/2*(v*(2*cm_k/h)/((2*cm_k/h) + cm_h) + k*difer -v * 1/2 - k/h);
                        A(j,j+1) = -v * 1/2 + k/h;
                    end
                end
            case cant_celdas
                %Caso ultima celda
                if (cbd_d ~= -1) %Condicion Dirichlet
                    A(j,j-1) = 1/2 * (v * 1/2 + k/h);
                    A(j,j) = h/dt + 1/2 *(v * 1/2 + k/h * -3);
                else
                    if(cbn_d ~= -1)
                        %Condicion Neumann
                        A(j,j-1) = 1/2*(v * 1/2 + k/h);
                        A(j,j) = h/dt + 1/2*(-v + v * 1/2 - k/h);
                    else
                        %Condicion Mixta
                        difer2= 2/h - 4*cm_k/(2*h*cm_k+cm_h*h^2);
                        A(j,j) = h/dt + 1/2*(v*(2*cm_k/h)/((-2*cm_k/h) + cm_h) + k*difer2 + v * 1/2 - k/h);
                        A(j,j+1) = 1/2*(v * 1/2 + k/h);
                    end
                end
            otherwise
                %celdas interna
                A(j,j-1) = 1/2*(v * 1/2 + k/h);
                A(j,j) = h/dt + 1/2 *(k/h * -2);
                A(j,j+1) = 1/2*(-v * 1/2 + k/h);
        end
    end
    
    
    %vector del lado izquierdo (todos los terminos que no tienen icognita fi y condiciones de borde)
    b = zeros(cant_celdas,1);
    
    for j = 1:cant_celdas
        %Partes de la ecuacion de transmision que colaboran con la matriz
        switch j
            case 1
                %Caso que es la primera celda
                if (cbd_i ~= -1)
                    %Condicion Dirichlet
                    ant1= -1/2*(- v * cbd_i  - (2*k/h)*(cbd_i -temp_t(j,i)) + (v/2)*(temp_t(j,i)+temp_t(j+1,i)) - (k/h)*(temp_t(j+1,i)-temp_t(j,i)));
                    b(j) = -Q*h - k/h * cbd_i - (1/2)*v * cbd_i + temp_t(j,i)*h/dt + ant1;%Parte de la actual mas parte del tiempo anterior
                else
                    if (cbn_i ~= -1)
                        %Condicion Neumann
                        ant2= -1/2*( -v*(cbn_i*(-h/2)+temp_t(j,i))  - k*cbn_i + (v/2)*(temp_t(j,i)+temp_t(j+1,i)) - (k/h)*(temp_t(j+1,i)-temp_t(j,i)));
                        b(j) = -Q*h  - (1/2)*(-k*cbn_i - v*(-1*(h/2)*cbn_i))+ temp_t(j,i)*h/dt + ant2;%el termino advectivo actua en la direccion de la cara
                    else
                        %Condicion Mixta
                        difer= 2*cm_h*cm_finf/(2*cm_k+cm_h*h) + (4*cm_k*temp(j,i))/(2*h*cm_k+cm_h*h^2) - 2*temp(j,i)/h;
                        ant3= -1/2*(-v*(cm_h*cm_finf+(2*cm_k/h)*temp(j,i))/((2*cm_k/h)+cm_h)  - k*difer  - k*cbn_i + (v/2)*(temp_t(j,i)+temp_t(j+1,i)) - (k/h)*(temp_t(j+1,i)-temp_t(j,i)));
                        b(j) = -Q*h - (1/2)*(-v*cm_finf*(cm_h)/((2*cm_k/h)+cm_h) - k*2*cm_h*cm_finf/(2*cm_k+cm_h*h)) + temp_t(j,i)*h/dt + ant3;
                    end
                end
            case cant_celdas
                %Caso ultima celda
                if (cbd_d ~= -1)
                    %Condicion Dirichlet
                    ant1= -1/2*(-(v/2)*(temp_t(j-1,i)+temp_t(j,i)) - (k/h)*(temp_t(j-1,i)-temp_t(j,i)) + v * cbd_d  - (2*k/h)*(cbd_d -temp_t(j,i)));
                    b(j) = -Q*h  - k/h * cbd_d + (1/2)*v*cbd_d + temp_t(j,i)*h/dt + ant1;
                else
                    if (cbn_d ~= -1)
                        %Condicion Neumann
                        ant2= -1/2*(-(v/2)*(temp_t(j-1,i)+temp_t(j,i)) - (k/h)*(temp_t(j-1,i)-temp_t(j,i)) + v*(cbn_d*(h/2)+temp_t(j,i)) - k*cbn_d );
                        b(j) = -Q*h - (1/2)*(-k*cbn_d + v*(1*(h/2)*cbn_d)) + temp_t(j,i)*h/dt + ant2;%el termino advectivo actua en la direccion de la cara
                    else
                        %Condicion Mixta
                        difer2= 2*temp_t(j,i)/h - (4*cm_k*temp_t(j,i))/(2*h*cm_k+cm_h*h^2) - (2*cm_k*cm_finf)/(2*cm_k+cm_h*h);
                        ant3= -1/2*(-v*(cm_h*cm_finf+(2*cm_k/h)*temp(j,i))/((2*cm_k/h)+cm_h) - k*difer2  - k*cbn_i + (v/2)*(temp_t(j,i)+temp_t(j-1,i)) - (k/h)*(temp_t(j-1,i)-temp_t(j,i)));
                        b(j) = -Q*h  (1/2)*(v*cm_finf*(cm_h)/((-2*cm_k/h)+cm_h) - k*(-2*cm_h*cm_finf/(2*cm_k+cm_h*h))) + temp_t(j,i)*h/dt + ant3;
                    end
                end
            otherwise
                %celdas interna
                ant4= -1/2*(-(v/2)*(temp_t(j-1,i)+temp_t(j,i)) - (k/h)*(temp_t(j-1,i)-temp_t(j,i)) + (v/2)*(temp_t(j+1,i)+temp_t(j,i)) - (k/h)*(temp_t(j+1,i)-temp_t(j,i)));
                b(j) = -Q*h + temp_t(i)*h/dt + ant4;
        end
    end
    %Resolucion del sistema
    temp=A\b;
    temp_t = [temp_t(:,:), temp(:)];
end

%plot(temp_t);
