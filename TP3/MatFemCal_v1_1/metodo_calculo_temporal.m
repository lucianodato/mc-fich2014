function [u] = metodo_calculo_temporal(K, C, f, dt, u_anterior,flag_metodo)
switch flag_metodo
    case 0 %forward
        u = C\(dt*(( C*1/dt + K)*u_anterior - f));
    case 1 %backward
        KC = K + 1/dt * C;
        u = KC\( f + 1/dt * C * u_anterior);
    otherwise %CN
        KC = K/2 - 1/dt * C;
        u= KC\(f + - (K*u_anterior)/2 - (C*u_anterior)/dt); % f = (fn_1 + fn_0)/2
end
