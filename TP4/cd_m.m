function [ temp ] = cd(cant_celdas,ini,fin,k,Q,v,cm_h,cm_k,cm_finf,cbd_i,cbd_d,cbn_i,cbn_d)
%Resolucion estacionaria por metodo de volumen finito para la ecuacion de
%transmision de Calor en 1D con termino advectivo por central difference y
%con celdas de tamaño regular

L = fin - ini;
h = L/cant_celdas;%distancia que separa las celdas (si son regulares)

%Matriz del lado derecho (todos los terminos que tienen incognitas fi)
A = zeros(cant_celdas,cant_celdas);

for i = 1:cant_celdas
    %Partes de la ecuacion de transmision que colaboran con la matriz
    switch i
        case 1
            %Caso que es la primera celda
            if (cbd_i ~= -1) %Condicion Dirichlet
                A(i,i) = -v * 1/2 + k/h * -3;
                A(i,i+1) = -v * 1/2 + k/h;
            else
                if(cbn_i ~= -1)
                    %Condicion Neumann
                    A(i,i) = v - v * 1/2 + k/h;
                    A(i,i+1) = -v * 1/2 + k/h;
                else
                    %Condicion Mixta
                    difer=4*cm_k/(2*h*cm_k+cm_h*h^2) - 2/h;
                    A(i,i) = v*(2*cm_k/h)/((2*cm_k/h) + cm_h) + k*difer + -v * 1/2 - k/h;
                    A(i,i+1) = -v * 1/2 + k/h;
                end
            end
        case cant_celdas
            %Caso ultima celda
            if (cbd_d ~= -1) %Condicion Dirichlet
                A(i,i-1) = v * 1/2 + k/h;
                A(i,i) = v * 1/2 + k/h * -3;
            else
                if(cbn_d ~= -1)
                    %Condicion Neumann
                    A(i,i-1) = v * 1/2 + k/h;
                    A(i,i) = -v + v * 1/2 - k/h ;
                else
                    %Condicion Mixta
                    difer2= 2/h - 4*cm_k/(2*h*cm_k+cm_h*h^2);
                    A(i,i) = v*(2*cm_k/h)/((2*cm_k/h) + cm_h) + k*difer2 + v * 1/2 - k/h;
                    A(i,i-1) = v * 1/2 - k/h;
                end
            end
        otherwise
            %celdas interna
            A(i,i-1) = v * 1/2 + k/h;
            A(i,i) = k/h * -2;
            A(i,i+1) = -v * 1/2 + k/h;
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
                b(i) = -Q*h - 2*k/h * cbd_i - v * cbd_i;
            else
                if (cbn_i ~= -1)
                    %Condicion Neumann
                    b(i) = -Q*h - k*cbn_i - v * (-1*(h/2)*cbn_i);%el termino advectivo actua en la direccion de la cara
                else
                    %Condicion Mixta
                    b(i) = -Q*h - v*cm_finf*(cm_h)/((2*cm_k/h)+cm_h) - k*2*cm_h*cm_finf/(2*cm_k+cm_h*h);
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
                    b(i) = -Q*h - k*cbn_d + v * (1*(h/2)*cbn_d);%el termino advectivo actua en la direccion de la cara
                else
                    %Condicion Mixta
                    b(i) = -Q*h + v*cm_finf*(cm_h)/((2*cm_k/h)+cm_h) - k*(-2*cm_h*cm_finf/(2*cm_k+cm_h*h)) ;
                end
            end
        otherwise
            %celdas interna
            b(i) = -Q*h;
    end
end

%Resolucion del sistema
temp=A\b;

end

