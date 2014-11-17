%codigo para el calculo de tensiones especificadas en un punto a partir de
%un sistema ya calculado y con valores en los nodos conocidos
syms x y;
punto_calculo = [0.75,0.75];

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
for i=1:length(elementos_cuadrangulares)-1
    x_max=max(desplazamientos(elementos_cuadrangulares(i,:),1));
    x_min=min(desplazamientos(elementos_cuadrangulares(i,:),1));
    y_max=max(desplazamientos(elementos_cuadrangulares(i,:),2));
    y_min=min(desplazamientos(elementos_cuadrangulares(i,:),2));
    %Si esta dentro del elemento entonces es el elemento para interpolar
    if(punto_calculo(1)<x_max && punto_calculo(1)>x_min && punto_calculo(2) <y_max && punto_calculo(2) >y_min)
        elemento_c = i;
    else
        elemento_c = -1;%por si no esta dentro de ningun cuadrangulo
    end
end

%Reocorro los trianguares
for i=1:length(elementos_triangulares)-1
    x_max=max(desplazamientos(elementos_triangulares(i,:),1));
    x_min=min(desplazamientos(elementos_triangulares(i,:),1));
    y_max=max(desplazamientos(elementos_triangulares(i,:),2));
    y_min=min(desplazamientos(elementos_triangulares(i,:),2));
    %Si esta dentro del elemento entonces es el elemento para interpolar
    if(punto_calculo(1)<x_max && punto_calculo(1)>x_min && punto_calculo(2) <y_max && punto_calculo(2) >y_min)
        elemento_t = i;
    else
        elemento_t = -1;%por si no esta dentro de ningun trirangulo
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
    
    l1= 1/2*sqrt((xj-xi)^2 + (yj-yi)^2);
    l2= 1/2*sqrt((xk-xj)^2 + (yk-yj)^2);
    A = l1*l2;
    
    %Familia de forma cuadrangular
    N1 = ((l1-x)*(l2-y))/(4*l1*l2);
    N2 = ((l1+x)*(l2-y))/(4*l1*l2);
    N3 = ((l1+x)*(l2+y))/(4*l1*l2);
    N4 = ((l1-x)*(l2+y))/(4*l1*l2);
    
    u1p=desplazamientos(elementos_cuadrangulares(elemento_c,1),3);
    u2p=desplazamientos(elementos_cuadrangulares(elemento_c,2),3);
    u3p=desplazamientos(elementos_cuadrangulares(elemento_c,3),3);
    u4p=desplazamientos(elementos_cuadrangulares(elemento_c,4),3);
    v1p=desplazamientos(elementos_cuadrangulares(elemento_c,1),4);
    v2p=desplazamientos(elementos_cuadrangulares(elemento_c,2),4);
    v3p=desplazamientos(elementos_cuadrangulares(elemento_c,3),4);
    v4p=desplazamientos(elementos_cuadrangulares(elemento_c,4),4);
    
    N1p = subs(N1,[x y],[punto_calculo(1),punto_calculo(2)]);
    N2p = subs(N2,[x y],[punto_calculo(1),punto_calculo(2)]);
    N3p = subs(N3,[x y],[punto_calculo(1),punto_calculo(2)]);
    N4p = subs(N4,[x y],[punto_calculo(1),punto_calculo(2)]);
    
    if(N1p+N2p+N3p+N4p == 1)
        valor_u = u1p*N1p + u2p*N2p + u3p*N3p + u4p*N4p;
        valor_v = v1p*N1p + v2p*N2p + v3p*N3p + v4p*N4p;
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
    else
        disp('Algo esta mal, las funciones de forma no suman 1');
    end
end
