%codigo para el calculo de temperaturas especificadas en un punto a partir de
%un sistema ya calculado y con valores en los nodos conocidos
syms x y;
punto_calculo = [1.75,0.75];

%Variables del problema
k=1;

%Tabla con valores cargados
% coordenada x - coordenada y - temperatura en el nodo
temperaturas = [
    0.0, 0.0, 0.0;
    1.0, 0.0, 0.7189;
    2.0, 0.0, -0.05839;
    2.0, 0.5, 0.55891;
    1.25, 0.5, 0.40925;
    0.0, 0.5, 0.0;
    0.0, 1.0, 0.0;
    1.0, 1.0, -0.13365;
    2.0, 1.0, 0.17009];

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
    x_max=max(temperaturas(elementos_cuadrangulares(i,:),1));
    x_min=min(temperaturas(elementos_cuadrangulares(i,:),1));
    y_max=max(temperaturas(elementos_cuadrangulares(i,:),2));
    y_min=min(temperaturas(elementos_cuadrangulares(i,:),2));
    %Si esta dentro del elemento entonces es el elemento para interpolar
    if(punto_calculo(1)<=x_max && punto_calculo(1)>=x_min && punto_calculo(2) <=y_max && punto_calculo(2) >=y_min)
        elemento_c = i;
    end
end

%Reocorro los trianguares
elemento_t = -1;%por si no esta dentro de ningun trirangulo
for i=1:size(elementos_triangulares,1)
    xi=temperaturas(elementos_triangulares(i,1),1);
    yi=temperaturas(elementos_triangulares(i,1),2);
    xj=temperaturas(elementos_triangulares(i,2),1);
    yj=temperaturas(elementos_triangulares(i,2),2);
    xk=temperaturas(elementos_triangulares(i,3),1);
    yk=temperaturas(elementos_triangulares(i,3),2);
    
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
    xi=temperaturas(elementos_cuadrangulares(elemento_c,1),1);
    yi=temperaturas(elementos_cuadrangulares(elemento_c,1),2);
    xj=temperaturas(elementos_cuadrangulares(elemento_c,2),1);
    yj=temperaturas(elementos_cuadrangulares(elemento_c,2),2);
    xk=temperaturas(elementos_cuadrangulares(elemento_c,3),1);
    yk=temperaturas(elementos_cuadrangulares(elemento_c,3),2);
    xl=temperaturas(elementos_cuadrangulares(elemento_c,4),1);
    yl=temperaturas(elementos_cuadrangulares(elemento_c,4),2);
    
    t1p=temperaturas(elementos_cuadrangulares(elemento_c,1),3);
    t2p=temperaturas(elementos_cuadrangulares(elemento_c,2),3);
    t3p=temperaturas(elementos_cuadrangulares(elemento_c,3),3);
    t4p=temperaturas(elementos_cuadrangulares(elemento_c,4),3);
    
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
    dN1x=diff(N1,x);
    dN1y=diff(N1,y);
    dN2x=diff(N2,x);
    dN2y=diff(N2,y);
    dN3x=diff(N3,x);
    dN3y=diff(N3,y);
    dN4x=diff(N4,x);
    dN4y=diff(N4,y);
        
    N1p = subs(N1,[x y],[punto_calculo(1),punto_calculo(2)]);
    N2p = subs(N2,[x y],[punto_calculo(1),punto_calculo(2)]);
    N3p = subs(N3,[x y],[punto_calculo(1),punto_calculo(2)]);
    N4p = subs(N4,[x y],[punto_calculo(1),punto_calculo(2)]);
    dN1xp= subs(dN1x,[x y],[punto_calculo(1),punto_calculo(2)]);
    dN1yp= subs(dN1y,[x y],[punto_calculo(1),punto_calculo(2)]);
    dN2xp= subs(dN2x,[x y],[punto_calculo(1),punto_calculo(2)]);
    dN2yp= subs(dN2y,[x y],[punto_calculo(1),punto_calculo(2)]);
    dN3xp= subs(dN3x,[x y],[punto_calculo(1),punto_calculo(2)]);
    dN3yp= subs(dN3y,[x y],[punto_calculo(1),punto_calculo(2)]);
    dN4xp= subs(dN4x,[x y],[punto_calculo(1),punto_calculo(2)]);
    dN4yp= subs(dN4y,[x y],[punto_calculo(1),punto_calculo(2)]);
    
    if(N1p+N2p+N3p+N4p == 1)
        valor_temperatura = t1p*N1p + t2p*N2p + t3p*N3p + t4p*N4p;
        valor_gradiente = -k*[t1p*dN1xp + t2p*dN2xp + t3p*dN3xp + t4p*dN4xp,t1p*dN1yp + t2p*dN2yp + t3p*dN3yp + t4p*dN4yp];
    else
        disp('Algo esta mal, las funciones de forma no suman 1');
    end
end



%-----------------Si es un triangulo--------------------
if(elemento_t ~=-1)
    xi=temperaturas(elementos_triangulares(elemento_t,1),1);
    yi=temperaturas(elementos_triangulares(elemento_t,1),2);
    xj=temperaturas(elementos_triangulares(elemento_t,2),1);
    yj=temperaturas(elementos_triangulares(elemento_t,2),2);
    xk=temperaturas(elementos_triangulares(elemento_t,3),1);
    yk=temperaturas(elementos_triangulares(elemento_t,3),2);
    
    Area = det([1 xi yi; 1 xj yj; 1 xk yk])/2;
    
    ai=(xj*yk-xk*yj); bi=(yj-yk); ci=(xk-xj);
    aj=(xk*yi-xi*yk); bj=(yk-yi); cj=(xi-xk);
    ak=(xi*yj-xj*yi); bk=(yi-yj); ck=(xj-xi);
    
    %Familia de forma triangular
    N1=(1/(2*Area))*(ai + bi*x + ci*y);
    N2=(1/(2*Area))*(aj + bj*x + cj*y);
    N3=(1/(2*Area))*(ak + bk*x + ck*y);
    dN1x=diff(N1,x);
    dN1y=diff(N1,y);
    dN2x=diff(N2,x);
    dN2y=diff(N2,y);
    dN3x=diff(N3,x);
    dN3y=diff(N3,y);
    
    t1p=temperaturas(elementos_triangulares(elemento_t,1),3);
    t2p=temperaturas(elementos_triangulares(elemento_t,2),3);
    t3p=temperaturas(elementos_triangulares(elemento_t,3),3);
    
    N1p = subs(N1,[x y],[punto_calculo(1),punto_calculo(2)]);
    N2p = subs(N2,[x y],[punto_calculo(1),punto_calculo(2)]);
    N3p = subs(N3,[x y],[punto_calculo(1),punto_calculo(2)]);
    dN1xp= subs(dN1x,[x y],[punto_calculo(1),punto_calculo(2)]);
    dN1yp= subs(dN1y,[x y],[punto_calculo(1),punto_calculo(2)]);
    dN2xp= subs(dN2x,[x y],[punto_calculo(1),punto_calculo(2)]);
    dN2yp= subs(dN2y,[x y],[punto_calculo(1),punto_calculo(2)]);
    dN3xp= subs(dN3x,[x y],[punto_calculo(1),punto_calculo(2)]);
    dN3yp= subs(dN3y,[x y],[punto_calculo(1),punto_calculo(2)]);
    
    if(N1p+N2p+N3p == 1)
        valor_temperatura = t1p*N1p + t2p*N2p + t3p*N3p;
        valor_gradiente = -k*[t1p*dN1xp + t2p*dN2xp + t3p*dN3xp,t1p*dN1yp + t2p*dN2yp + t3p*dN3yp];
    else
        disp('Algo esta mal, las funciones de forma no suman 1');
    end
end
