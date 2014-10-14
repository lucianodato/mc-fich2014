% Las entradas vienen en este formato
%
% coord =
% 
%    40.0000   54.2857
%    40.0000   51.4286
%    42.5981   52.8571
% 
% 
% dmat =
% 
%     63     0
%      0    63
% 
% 
% heat =
% 
%      0

function Elem_C_Mat = obtener_C_mat_quad(nodes, rho, cp)

fform = @(s,t)[(1-s-t+s*t)/4,(1+s-t-s*t)/4,(1+s+t+s*t)/4,(1-s+t-s*t)/4];

pospg = [ -0.577350269189626E+00 , 0.577350269189626E+00 ];
pespg = [  1.0E+00 , 1.0E+00]; % Esta variable es el peso de wi y wj, pesos relativos a nu y chi en la integración
M = zeros(4,4);

for i=1 : 2
    for j=1 : 2
        lcffm = fform(pospg(i),pospg(j)) ;    % FF at gauss point
        xjacm = lcffm*nodes ;                 % Jacobian matrix
        ctder = xjacm\lcffm ;                 % FF Cartesian derivates
        darea = det(xjacm)*pespg(i)*pespg(j);
        
        bmat = [];
        for inode = 1 : 4
            bmat = [ bmat , [ctder(1,inode);
                ctder(2,inode) ] ] ;
        end
        
        M = M + (transpose(bmat)*bmat)*darea;
        
    end
end

Elem_C_Mat = rho * cp * M;

end