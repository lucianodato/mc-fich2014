%Ejercicio 1 
%-----------DATOS-------------
paso = 1/4;%que tan fina es la subdivision para los caluculos
b1 = 0;%punto inicial o borde 1
b2 = 1;%punto final o borde 2
inciso = 3; %(1->a 2->b 3->c)
%-----------------------------

p = b1:paso:b2;
L = max(p);%Maximo punto para la normalizacion de furier
cant_puntos = length(p);%Cantidad total de puntos segun b1 b2 y paso

syms x;

Klm = zeros(cant_puntos,cant_puntos);

switch inciso
    case 1
        for i = 1:cant_puntos
            for j = 1:cant_puntos
                %Calculamos Ni y Nj
                if(i==j || i==j-1 || i==j+1)
                    Ni = 1;
                    Nj = 1;
                else
                    Ni = 0;
                    Nj = 0;
                end
                %Calculamos Klm
                Klm(i,j) = int(Ni*Nj,x,0,1);
            end
        end
    case 2
        for i = 1:cant_puntos
            for j = 1:cant_puntos
                %Calculamos Ni y Nj
                Ni = x^i;
                Nj = x^j;
                %Calculamos Klm
                Klm(i,j) = int(Ni*Nj,x,0,1);
            end
        end
    case 3
        l = [1,5,3,4,2];
        for i = 1:cant_puntos
            for j = 1:cant_puntos
                %Calculamos Ni y Nj
                if(i==j || i==j-1 || i==j+1)
                    Ni = 1;
                    Nj = 1;
                else
                    Ni = 0;
                    Nj = 0;
                end
                %Calculamos Klm
                Klm(l(i),l(j)) = int(Ni*Nj,x,0,1);
            end
        end
end

Klm
