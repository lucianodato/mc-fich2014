%codigo para el calculo de tensiones especificadas en un punto a partir de
%un sistema ya calculado y con valores en los nodos conocidos
syms x y;
punto_calculo = [1.5,0.75];

%Variables del problema
E=2.1e9;%Modulo de Young
v=0.3;%Poisson
t=0.05;%espesor
%Matriz D tension plana
Dt  = E/(1-v*v)*[1, v, 0;
                  v, 1, 0;
                  0, 0, 0.5*(1-v)];

%Tabla con valores cargados
% coordenada x - coordenada y - desplazamiento u - desplazamiento v
desplazamientos = [
    0.0, 0.0, 0.0, 0.0;
    1.0, 0.0, 0.7189, 0.66431;
    2.0, 0.0, -0.05839, 0.0;
    2.0, 0.5, 0.55891, 1.7722;
    1.25, 0.5, 0.40925, 0.89222;
    0.0, 0.5, 0.0, 0.0;
    0.0, 1.0, 0.0, 0.0;
    1.0, 1.0, -0.13365, 0.45779;
    2.0, 1.0, 0.17009, 3.9252];

elementos_cuadrangulares = [
    1,2,5,6;
    2,3,4,5;
    5,4,9,8];

elementos_triangulares = [
    6,5,8;
    6,8,7];

%En primera instancia se determina en que elemento hay que calcular el
%resultado

%Reocorro los cuadranguares
elemento_c = -1;%por si no esta dentro de ningun cuadrangulo
for i=1:size(elementos_cuadrangulares,1)
    x_max=max(desplazamientos(elementos_cuadrangulares(i,:),1));
    x_min=min(desplazamientos(elementos_cuadrangulares(i,:),1));
    y_max=max(desplazamientos(elementos_cuadrangulares(i,:),2));
    y_min=min(desplazamientos(elementos_cuadrangulares(i,:),2));
    %Si esta dentro del elemento entonces es el elemento para interpolar
    if(punto_calculo(1)<=x_max && punto_calculo(1)>=x_min && punto_calculo(2)<=y_max && punto_calculo(2)>=y_min)
        elemento_c = i;
    end
end

%Reocorro los trianguares
elemento_t = -1;%por si no esta dentro de ningun trirangulo
for i=1:size(elementos_triangulares,1)
    xi=desplazamientos(elementos_triangulares(i,1),1);
    yi=desplazamientos(elementos_triangulares(i,1),2);
    xj=desplazamientos(elementos_triangulares(i,2),1);
    yj=desplazamientos(elementos_triangulares(i,2),2);
    xk=desplazamientos(elementos_triangulares(i,3),1);
    yk=desplazamientos(elementos_triangulares(i,3),2);
    
    %calculo las orientaciones para determinar si esta dentro
    orient=sign((xi-xk)*(yj-yk)-(yi-yk)*(xj-xk));
    orient1=sign((xi-punto_calculo(1))*(yj-punto_calculo(2))-(yi-punto_calculo(2))*(xj-punto_calculo(1)));
    orient2=sign((xj-punto_calculo(1))*(yk-punto_calculo(2))-(yj-punto_calculo(2))*(xk-punto_calculo(1)));
    orient3=sign((xk-punto_calculo(1))*(yi-punto_calculo(2))-(yk-punto_calculo(2))*(xi-punto_calculo(1)));
    %si las orientaciones son iguales entonces el punto esta dentro del
    %triangulo
    if(orient>0)%orientacion positiva
        if((orient1+orient2+orient3)>=3)
            elemento_t=i;
        end
    end
    if(orient<0)%orientacion negativa
        if((orient1+orient2+orient3)<=-3)
            elemento_t=i;
        end
    end
end

%-----------------Si es un cuadrangulo--------------------
if(elemento_c ~=-1)
    xi=desplazamientos(elementos_cuadrangulares(elemento_c,1),1);
    yi=desplazamientos(elementos_cuadrangulares(elemento_c,1),2);
    xj=desplazamientos(elementos_cuadrangulares(elemento_c,2),1);
    yj=desplazamientos(elementos_cuadrangulares(elemento_c,2),2);
    xk=desplazamientos(elementos_cuadrangulares(elemento_c,3),1);
    yk=desplazamientos(elementos_cuadrangulares(elemento_c,3),2);
    xl=desplazamientos(elementos_cuadrangulares(elemento_c,4),1);
    yl=desplazamientos(elementos_cuadrangulares(elemento_c,4),2);
        
    u1p=desplazamientos(elementos_cuadrangulares(elemento_c,1),3);
    u2p=desplazamientos(elementos_cuadrangulares(elemento_c,2),3);
    u3p=desplazamientos(elementos_cuadrangulares(elemento_c,3),3);
    u4p=desplazamientos(elementos_cuadrangulares(elemento_c,4),3);
    v1p=desplazamientos(elementos_cuadrangulares(elemento_c,1),4);
    v2p=desplazamientos(elementos_cuadrangulares(elemento_c,2),4);
    v3p=desplazamientos(elementos_cuadrangulares(elemento_c,3),4);
    v4p=desplazamientos(elementos_cuadrangulares(elemento_c,4),4);
    
    %Usando funciones de forma triangulares en 2d (piramides)
    A = [1 xi yi xi*yi;
        1 xj yj xj*yj;
        1 xk yk xk*yk;
        1 xl yl xl*yl];
    f1 = [1 0 0 0]';
    f2 = [0 1 0 0]';
    f3 = [0 0 1 0]';
    f4 = [0 0 0 1]';
    
    N1i=A\f1;
    N2i=A\f2;
    N3i=A\f3;
    N4i=A\f4;
    
    %Familia de forma cuadrangular
    N1 = N1i(1)+N1i(2)*x+N1i(3)*y+N1i(4)*x*y;
    N2 = N2i(1)+N2i(2)*x+N2i(3)*y+N2i(4)*x*y;
    N3 = N3i(1)+N3i(2)*x+N3i(3)*y+N3i(4)*x*y;
    N4 = N4i(1)+N4i(2)*x+N4i(3)*y+N4i(4)*x*y;
    
    B1=[N1i(2)+N1i(4)*y 0; 0 N1i(3)+N1i(4)*x;N1i(3)+N1i(4)*x N1i(2)+N1i(4)*y];
    B2=[N2i(2)+N2i(4)*y 0; 0 N2i(3)+N2i(4)*x;N2i(3)+N2i(4)*x N2i(2)+N2i(4)*y];
    B3=[N3i(2)+N3i(4)*y 0; 0 N3i(3)+N3i(4)*x;N3i(3)+N3i(4)*x N3i(2)+N3i(4)*y];
    B4=[N4i(2)+N4i(4)*y 0; 0 N4i(3)+N4i(4)*x;N4i(3)+N4i(4)*x N4i(2)+N4i(4)*y];
        
    N1p = subs(N1,[x y],[punto_calculo(1),punto_calculo(2)]);
    N2p = subs(N2,[x y],[punto_calculo(1),punto_calculo(2)]);
    N3p = subs(N3,[x y],[punto_calculo(1),punto_calculo(2)]);
    N4p = subs(N4,[x y],[punto_calculo(1),punto_calculo(2)]);
    B1p = subs(B1,[x y],[punto_calculo(1),punto_calculo(2)]);
    B2p = subs(B2,[x y],[punto_calculo(1),punto_calculo(2)]);
    B3p = subs(B3,[x y],[punto_calculo(1),punto_calculo(2)]);
    B4p = subs(B4,[x y],[punto_calculo(1),punto_calculo(2)]);
    
    if(N1p+N2p+N3p+N4p == 1)
        valor_u = u1p*N1p + u2p*N2p + u3p*N3p + u4p*N4p;
        valor_v = v1p*N1p + v2p*N2p + v3p*N3p + v4p*N4p;
        deformacion = B1p*[u1p;v1p] + B2p*[u2p;v2p] + B3p*[u3p;v3p] + B4p*[u4p;v4p];
        tension=Dt*deformacion;
    else
        disp('Algo esta mal, las funciones de forma no suman 1');
    end
end

%-----------------Si es un triangulo--------------------
if(elemento_t ~=-1)
    xi=desplazamientos(elementos_triangulares(elemento_t,1),1);
    yi=desplazamientos(elementos_triangulares(elemento_t,1),2);
    xj=desplazamientos(elementos_triangulares(elemento_t,2),1);
    yj=desplazamientos(elementos_triangulares(elemento_t,2),2);
    xk=desplazamientos(elementos_triangulares(elemento_t,3),1);
    yk=desplazamientos(elementos_triangulares(elemento_t,3),2);
    
    Area = det([1 xi yi; 1 xj yj; 1 xk yk])/2;
    
    ai=(xj*yk-xk*yj); bi=(yj-yk); ci=(xk-xj);
    aj=(xk*yi-xi*yk); bj=(yk-yi); cj=(xi-xk);
    ak=(xi*yj-xj*yi); bk=(yi-yj); ck=(xj-xi);
    
    %Familia de forma triangular
    N1=(1/(2*Area))*(ai + bi*x + ci*y);
    N2=(1/(2*Area))*(aj + bj*x + cj*y);
    N3=(1/(2*Area))*(ak + bk*x + ck*y);
    B1=(1/(2*Area))*[bi 0; 0 ci; ci bi];
    B2=(1/(2*Area))*[bj 0; 0 cj; cj bj];
    B3=(1/(2*Area))*[bk 0; 0 ck; ck bk];
        
    u1p=desplazamientos(elementos_triangulares(elemento_t,1),3);
    u2p=desplazamientos(elementos_triangulares(elemento_t,2),3);
    u3p=desplazamientos(elementos_triangulares(elemento_t,3),3);
    v1p=desplazamientos(elementos_triangulares(elemento_t,1),4);
    v2p=desplazamientos(elementos_triangulares(elemento_t,2),4);
    v3p=desplazamientos(elementos_triangulares(elemento_t,3),4);
    
    N1p = subs(N1,[x y],[punto_calculo(1),punto_calculo(2)]);
    N2p = subs(N2,[x y],[punto_calculo(1),punto_calculo(2)]);
    N3p = subs(N3,[x y],[punto_calculo(1),punto_calculo(2)]);
    
    if(N1p+N2p+N3p == 1)
        valor_u = u1p*N1p + u2p*N2p + u3p*N3p;
        valor_v = v1p*N1p + v2p*N2p + v3p*N3p;
        deformacion = B1*[u1p;v1p] + B2*[u2p;v2p] + B3*[u3p;v3p];
        tension=Dt*deformacion;
    else
        disp('Algo esta mal, las funciones de forma no suman 1');
    end
end
