 %ejer 1
 x = [0:5];
 l=1;
 %estos a y b se hallaron a partir del metodo para resolver la edo y las condiciones de borde
 a = -exp(-l)/(exp(-l) - exp(l));
 b = exp(l)/(exp(-l) - exp(l));
 t = a.*exp(x) + b.*exp(-x) + 1;
 plot(t);