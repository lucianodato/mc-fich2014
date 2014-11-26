cant = 5;
t_max = 0.01;
t_ini=0;
dt=0.002;
phi0=0;
dif_phi6 = 1;
cant_pasos_tiempo = (t_max-t_ini)/dt;
b= zeros(cant,1);

A= [132 27.5 0 0 0;
    -27.5 105 22.5 0 0;
    0 -27.5 105 22.5 0;
    0 0 -27.5 105 22.5;
    0 0 0 -27.5 127.5];

%Phi tiempo 0.
phi1= 0;
phi2= 0;
phi3= 0;
phi4= 0;
phi5= 0;


temp_t = [phi1,phi2,phi3,phi4,phi5]';

for i=1:cant_pasos_tiempo
    
    st1=100*temp_t(1,i) - 0.5*(-100*phi0 - 10*(phi0 - temp_t(1,i)) + 50*(temp_t(1,i) + temp_t(2,i)) - 5*(temp_t(2,i) - temp_t(1,i)));
    st2=100*temp_t(2,i) - 0.5*(-50*(temp_t(1,i)+temp_t(2,i)) - 5*(temp_t(1,i) - temp_t(2,i)) + 50*(temp_t(2,i) + temp_t(3,i)) - 5*(temp_t(3,i) - temp_t(2,i)));
    st3=100*temp_t(3,i) - 0.5*(-50*(temp_t(2,i)+temp_t(3,i)) - 5*(temp_t(2,i) - temp_t(3,i)) + 50*(temp_t(3,i) + temp_t(4,i)) - 5*(temp_t(4,i) - temp_t(3,i)));
    st4=100*temp_t(4,i) - 0.5*(-50*(temp_t(3,i)+temp_t(4,i)) - 5*(temp_t(3,i) - temp_t(4,i)) + 50*(temp_t(4,i) + temp_t(5,i)) - 5*(temp_t(5,i) - temp_t(4,i)));
    st5=100*temp_t(5,i) - 0.5*(-50*(temp_t(4,i)+temp_t(5,i)) - 5*(temp_t(4,i) - temp_t(5,i)) + 10*dif_phi6 + 100*temp_t(5,i) -dif_phi6);
    
    b=[(1/5 + st1) (1/5 + st2) (1/5 + st3) (1/5 + st4) (-8.8 + st5)];
    
    a= A\b';
    temp_t = [temp_t(:,:),a(:)];
    
end

plot(temp_t);
legend('t0','t1','t2','t3','t4','t5');