%Definicion de parametros globales
cant_celdas = 5; %numeros de elementos en la recta 1D
ini = 0;%punto de inicio de la recta
fin = 1;%punto de fin de la recta
x=1:cant_celdas;

%Definicion de parametros especificos de la ecuacion de balance termico
k = -1;%Constante de difusividad
Q = -1;%fuente
v = -100;%velocidad
cm_h = 1;%h de la condicion mixta
cm_k = 1;%k de la condicion mixta si la hay
cm_finf = 1;%temperatura externa fi inf

%Definicion de las condiciones de borde (-1 significa que no aplica)
cbd_i = 0;%condicion de borde dirichlet izquierda
cbd_d = -1;%condicion de borde dirichlet derecha
cbn_i = -1;%condicion de borde neumann izquierda
cbn_d = 1;%condicion de borde neumann derecha

cdr = cd_m(cant_celdas,ini,fin,k,Q,v,cm_h,cm_k,cm_finf,cbd_i,cbd_d,cbn_i,cbn_d);
udr = ud_m(cant_celdas,ini,fin,k,Q,v,cm_h,cm_k,cm_finf,cbd_i,cbd_d,cbn_i,cbn_d);
anr = analitica(cant_celdas,ini,fin,-v,-Q);

plot(x,anr,'--*r',x,cdr,'--*g',x,udr,'--*b');
legend('Analitica','C.D.','U.D.');