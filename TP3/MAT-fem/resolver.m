cant = 3;
A= zeros(cant,cant);
b= zeros(cant,1);

A= [1 0;
    1/3 1/3];

b=[1 25/3];

a= A\b';

