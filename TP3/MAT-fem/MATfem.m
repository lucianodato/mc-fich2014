%% MAT-fem
%
% Clear memory and variables.
clear
% The variables are read as a MAT-fem subroutine
% pstrs = 1 indicate Plane Stress; 0 indicate Plane Strain
% young =   Young Modulus
% poiss =   Poission Ratio
% thick =   thickness only valid for Plane Stress
% denss =   density
% coordinates = [ x , y ] coordinate matrix nnode x ndime (2)
% elements    = [ inode, jnode, knode ] element connectivity  matrix
%               nelem x nnode; nnode = 3 for triangular elements and
%               nnode = 4 for quadrilateral elements
% fixnodes    = [node number, dimension, fixed value] matrix with
%               Dirichlet restrictions.
% pointload   = [node number, dimension, load value] matrix with
%               nodal loads.
% SideLoad    = [node number i, node number j, x load, y load] matrix with
%               line definition and uniform load applied in x and y
%               directions
% VariableSideLoad = [node number i, node number j,x load, y load,element] matrix 
%               with line definition and uniform load applied in x and y
%               directions
% midpointload   = [xpos , ypos, load x, load y,element] matrix with
%               nodal loads.

%Manual filling of variablesideload forces
syms x y;
x_function = (25000/2) * x - 10000.0; 
%y_function = y - 10.0; 
%variablesideload = [];
variablesideload = [
3 , 1 , 0.00000 , x_function, 10;
7 , 3 , 0.00000 , x_function, 12];

%Midpoint forces filling
midpointload = [];
%midpointload = [0.5 , 1.5 ,-1000.0 , 1000.0 , 5];
  
%Carga de archivo
file_name = input('Enter the file name :','s');

tic;                   % Start clock
ttim = 0;              % Initialize time counter

cd('Problems/');
eval (file_name);      % Read input file
cd('..');

% Finds basics dimensions
npnod  = size(coordinates,1);        % Number of nodes
nndof  = 2*npnod;                    % Number of total DOF
nelem  = size(elements,1);           % Number of elements
nnode  = size(elements,2);           % Number of nodes per element
neleq  = nnode*2;                    % Number of DOF per element

ttim = timing('Time needed to read the input file',ttim); %Reporting time

% Dimension the global matrices.
StifMat = sparse ( nndof , nndof );  % Create the global stiffness matrix
force   = sparse ( nndof , 1 );      % Create the global force vector

%  Material properties (Constant over the domain).
dmat = constt(young,poiss,pstrs);

ttim = timing('Time needed to  set initial values',ttim); %Reporting time

%  Element loop.
for ielem = 1 : nelem
    
    % Recover element properties
    lnods = elements(ielem,:);                        % Elem. connectivity
    coord(1:nnode,:) = coordinates(lnods(1:nnode),:); % Elem. coordinates
    
    % Evaluates the elemental stiffness matrix and mass force vector.
    if (nnode == 3)
        [ElemMat,ElemFor] = TrStif(coord,dmat ,thick,denss); % 3 Nds Triangle
    else
        [ElemMat,ElemFor] = QdStif(coord,dmat ,thick,denss); % 4 Nds Quad.
    end
    
    % Finds the equation number list for the i-th element
    eqnum = [];                                  % Clear the list
    for i =1 : nnode                             % Node cicle
        eqnum = [eqnum,lnods(i)*2-1,lnods(i)*2];   % Build the equation
    end                                          % number list
    
    % Assemble the force vector and the stiffness matrix
    for i = 1 : neleq
        force(eqnum(i)) = force(eqnum(i)) + ElemFor(i);
        for j = 1 : neleq
            StifMat(eqnum(i),eqnum(j)) = StifMat(eqnum(i),eqnum(j)) + ...
                ElemMat(i,j);
        end
    end
    
end  % End element loop

ttim = timing('Time to assamble the global system',ttim); %Reporting time

%  Add variabe side forces to the force vector (loads could be functions)
for i = 1 : size(variablesideload,1)
    syms x y;
    if (nnode == 3)
        %Area del triangulo (con las cordenadas de cada nodo de elemento)
        %element node
        node_1 = elements(variablesideload(i,5),1);
        node_2 = elements(variablesideload(i,5),2);
        node_3 = elements(variablesideload(i,5),3);
        
        xi = coordinates(node_1,1);
        yi = coordinates(node_1,2);
        xj = coordinates(node_2,1);
        yj = coordinates(node_2,2);
        xk = coordinates(node_3,1);
        yk = coordinates(node_3,2);
        
        %not used node
        if(variablesideload(i,1) ~= node_1 && variablesideload(i,2) ~= node_1)
            nu_node = 1;
        end
        if(variablesideload(i,1) ~= node_2 && variablesideload(i,2) ~= node_2)
            nu_node = 2;
        end
        if(variablesideload(i,1) ~= node_3 && variablesideload(i,2) ~= node_3)
            nu_node = 3;
        end

        A = 1/2* det([1 xi yi;1 xj yj;1 xk yk]);
        %Familia de forma triangular
        Ni = (1/2*A)* ((xj*yk-xk*yj) + (yj-yk)*x + (xk-xj)*y);%symbolic
        Nj = (1/2*A)* ((xk*yi-yk*xi) + (yk-yi)*x + (xi-xk)*y);
        Nk = (1/2*A)* ((xi*yj-yi*xj) + (yi-yj)*x + (xj-xi)*y);
        
        switch nu_node
            case 1 %node i is not used 
                ieqn = variablesideload(i,1)*2;         % Finds eq. number for the first node
                force(ieqn-1) = force(ieqn-1) + int(subs(Nj*variablesideload(i,3),x,xj),y,yj,yk);   % add x force
                force(ieqn  ) = force(ieqn  ) + int(subs(Nj*variablesideload(i,4),y,yj),x,xj,xk);   % add y force
                
                ieqn = variablesideload(i,2)*2;         % Finds eq. number for the second node
                force(ieqn-1) = force(ieqn-1) + int(subs(Nk*variablesideload(i,3),x,xk),y,yk,yj);   % add x force
                force(ieqn  ) = force(ieqn  ) + int(subs(Nk*variablesideload(i,4),y,yk),x,xk,xj);   % add y force
            case 2 %node j is not used 
                ieqn = variablesideload(i,1)*2;         % Finds eq. number for the first node
                force(ieqn-1) = force(ieqn-1) + int(subs(Ni*variablesideload(i,3),x,xi),y,yi,yk);   % add x force
                force(ieqn  ) = force(ieqn  ) + int(subs(Ni*variablesideload(i,4),y,yi),x,xi,xk);   % add y force
                
                ieqn = variablesideload(i,2)*2;         % Finds eq. number for the second node
                force(ieqn-1) = force(ieqn-1) + int(subs(Nk*variablesideload(i,3),x,xk),y,yk,yi);   % add x force
                force(ieqn  ) = force(ieqn  ) + int(subs(Nk*variablesideload(i,4),y,yk),x,xk,xi);   % add y force
            case 3 %node k is not used 
                ieqn = variablesideload(i,1)*2;         % Finds eq. number for the first node
                force(ieqn-1) = force(ieqn-1) + int(subs(Ni*variablesideload(i,3),x,xi),y,yi,yj);   % add x force
                force(ieqn  ) = force(ieqn  ) + int(subs(Ni*variablesideload(i,4),y,yi),x,xi,xj);   % add y force
                
                ieqn = variablesideload(i,2)*2;         % Finds eq. number for the second node
                force(ieqn-1) = force(ieqn-1) + int(subs(Nj*variablesideload(i,3),x,xj),y,yj,yi);   % add x force
                force(ieqn  ) = force(ieqn  ) + int(subs(Nj*variablesideload(i,4),y,yj),x,xj,xi);   % add y force
        end
    else
        %Area del cuadrangulo (con las cordenadas de cada nodo de elemento)
        
        %element node
        node_1 = elements(variablesideload(i,5),1);
        node_2 = elements(variablesideload(i,5),2);
        node_3 = elements(variablesideload(i,5),3);
        node_4 = elements(variablesideload(i,5),4);
        
        xi = coordinates(node_1,1);
        yi = coordinates(node_1,2);
        xj = coordinates(node_2,1);
        yj = coordinates(node_2,2);
        xk = coordinates(node_3,1);
        yk = coordinates(node_3,2);
        xl = coordinates(node_4,1);
        yl = coordinates(node_4,2);
        
        nu_node = [];
        %not used nodes
        if(variablesideload(i,1) ~= node_1 && variablesideload(i,2) ~= node_1 && variablesideload(i,1) ~= node_2 && variablesideload(i,2) ~= node_2)
            nu_node = [1 2];%i y j no usados
        end
        if(variablesideload(i,1) ~= node_1 && variablesideload(i,2) ~= node_1 && variablesideload(i,1) ~= node_4 && variablesideload(i,2) ~= node_4)
            nu_node = [1 4];%i y l no usados
        end
        if(variablesideload(i,1) ~= node_2 && variablesideload(i,2) ~= node_2 && variablesideload(i,1) ~= node_3 && variablesideload(i,2) ~= node_3)
            nu_node = [2 3];%j y k no usados
        end
        if(variablesideload(i,1) ~= node_3 && variablesideload(i,2) ~= node_3 && variablesideload(i,1) ~= node_4 && variablesideload(i,2) ~= node_4)
            nu_node = [3 4];%k y l no usados
        end
        
        %we have both not used nodes now we have to determine which side is
        %so what are those form function that do not contribute to the
        %force
        
        if (nu_node == [1 2])
            side = 3;%side 3 is used
        end
        if (nu_node == [2 3])
            side = 4;%side 4 is used
        end
        if (nu_node == [3 4])
            side = 1;%side 1 is used
        end
        if (nu_node == [1 4])
            side = 2;%side 2 is used
        end

        l1= 1/2*sqrt((xj-xi)^2 + (yj-yi)^2);
        l2= 1/2*sqrt((xk-xj)^2 + (yk-yj)^2);
        A = l1*l2;
        
        %Familia de forma cuadrangular
        Ni = ((l1-x)*(l2-y))/(4*l1*l2);%symbolic
        Nj = ((l1+x)*(l2-y))/(4*l1*l2);
        Nk = ((l1+x)*(l2+y))/(4*l1*l2);
        Nl = ((l1-x)*(l2+y))/(4*l1*l2);
        
        
        switch side
            case 1 %side 1 is used i and j are form functions
                ieqn = variablesideload(i,1)*2;         % Finds eq. number for the first node
                force(ieqn-1) = force(ieqn-1) + int(subs(Ni*variablesideload(i,3),x,xi),y,yi,yj);   % add x force
                force(ieqn  ) = force(ieqn  ) + int(subs(Ni*variablesideload(i,4),y,yi),x,xi,xj);   % add y force
                
                ieqn = variablesideload(i,2)*2;         % Finds eq. number for the second node
                force(ieqn-1) = force(ieqn-1) + int(subs(Nj*variablesideload(i,3),x,xj),y,yj,yi);   % add x force
                force(ieqn  ) = force(ieqn  ) + int(subs(Nj*variablesideload(i,4),y,yj),x,xj,xi);   % add y force
            case 2 %side 2 is used j and k are form functions
                ieqn = variablesideload(i,1)*2;         % Finds eq. number for the first node
                force(ieqn-1) = force(ieqn-1) + int(subs(Nj*variablesideload(i,3),x,xj),y,yj,yk);   % add x force
                force(ieqn  ) = force(ieqn  ) + int(subs(Nj*variablesideload(i,4),y,yj),x,xj,xk);   % add y force
                
                ieqn = variablesideload(i,2)*2;         % Finds eq. number for the second node
                force(ieqn-1) = force(ieqn-1) + int(subs(Nk*variablesideload(i,3),x,xk),y,yk,yj);   % add x force
                force(ieqn  ) = force(ieqn  ) + int(subs(Nk*variablesideload(i,4),y,yk),x,xk,xj);   % add y force
            case 3 %side 3 is used k and h are form functions
                ieqn = variablesideload(i,1)*2;         % Finds eq. number for the first node
                force(ieqn-1) = force(ieqn-1) + int(subs(Nk*variablesideload(i,3),x,xk),y,yk,yl);   % add x force
                force(ieqn  ) = force(ieqn  ) + int(subs(Nk*variablesideload(i,4),y,yk),x,xk,xl);   % add y force
                
                ieqn = variablesideload(i,2)*2;         % Finds eq. number for the second node
                force(ieqn-1) = force(ieqn-1) + int(subs(Nl*variablesideload(i,3),x,xl),y,yl,yk);   % add x force
                force(ieqn  ) = force(ieqn  ) + int(subs(Nl*variablesideload(i,4),y,yl),x,xl,xk);   % add y force
            case 4 %side 4 is used k and i are form functions
                ieqn = variablesideload(i,1)*2;         % Finds eq. number for the first node
                force(ieqn-1) = force(ieqn-1) + int(subs(Nl*variablesideload(i,3),x,xl),y,yl,yi);   % add x force
                force(ieqn  ) = force(ieqn  ) + int(subs(Nl*variablesideload(i,4),y,yl),x,xl,xi);   % add y force
                
                ieqn = variablesideload(i,2)*2;         % Finds eq. number for the second node
                force(ieqn-1) = force(ieqn-1) + int(subs(Ni*variablesideload(i,3),x,xi),y,yi,yl);   % add x force
                force(ieqn  ) = force(ieqn  ) + int(subs(Ni*variablesideload(i,4),y,yi),x,xi,xl);   % add y force
        end
    end
end

%  Add side forces to the force vector
for i = 1 : size(sideload,1)
    x=coordinates(sideload(i,1),:)-coordinates(sideload(i,2),:);
    l = sqrt(x*transpose(x));       % Finds the length of the side
    ieqn = sideload(i,1)*2;         % Finds eq. number for the first node
    force(ieqn-1) = force(ieqn-1) + l*sideload(i,3)/2;   % add x force
    force(ieqn  ) = force(ieqn  ) + l*sideload(i,4)/2;   % add y force
    
    ieqn = sideload(i,2)*2;         % Finds eq. number for the second node
    force(ieqn-1) = force(ieqn-1) + l*sideload(i,3)/2;   % add x force
    force(ieqn  ) = force(ieqn  ) + l*sideload(i,4)/2;   % add y force
end

%  Add midpoint loads conditions to the force vector
for i = 1 : size(midpointload,1)
    syms x y;
    if (nnode == 3)
        %Area del triangulo (con las cordenadas de cada nodo de elemento)
        xi = coordinates(elements(midpointload(i,5),1),1);
        yi = coordinates(elements(midpointload(i,5),1),2);
        xj = coordinates(elements(midpointload(i,5),2),1);
        yj = coordinates(elements(midpointload(i,5),2),2);
        xk = coordinates(elements(midpointload(i,5),3),1);
        yk = coordinates(elements(midpointload(i,5),3),2);
        
        A = 1/2* det([1 xi yi;1 xj yj;1 xk yk]);
        %Familia de forma triangular
        Ni = (1/2*A)* ((xj*yk-xk*yj) + (yj-yk)*x + (xk-xj)*y);
        Nj = (1/2*A)* ((xk*yi-yk*xi) + (yk-yi)*x + (xi-xk)*y);
        Nk = (1/2*A)* ((xi*yj-yi*xj) + (yi-yj)*x + (xj-xi)*y);
        N1= subs(Ni,[x,y],[midpointload(i,1),midpointload(i,2)]);
        N2= subs(Nj,[x,y],[midpointload(i,1),midpointload(i,2)]);
        N3= subs(Nk,[x,y],[midpointload(i,1),midpointload(i,2)]);
        
        %Para cada punto del elemento distribuyo la carga que se situa
        %entre medio del elemento
        ieqn = elements(midpointload(i,5),1)*2;         % Finds eq. number for the first node
        force(ieqn-1) = force(ieqn-1) + N1*midpointload(i,3)*A;   % add x force
        force(ieqn  ) = force(ieqn  ) + N1*midpointload(i,4)*A;   % add y force
        
        ieqn = elements(midpointload(i,5),2)*2;         % Finds eq. number for the second node
        force(ieqn-1) = force(ieqn-1) + N2*midpointload(i,3)*A;   % add x force
        force(ieqn  ) = force(ieqn  ) + N2*midpointload(i,4)*A;   % add y force
        
        ieqn = elements(midpointload(i,5),3)*2;         % Finds eq. number for the third node
        force(ieqn-1) = force(ieqn-1) + N3*midpointload(i,3)*A;   % add x force
        force(ieqn  ) = force(ieqn  ) + N3*midpointload(i,4)*A;   % add y force
        
    else
        %Area del cuadrangulo (con las cordenadas de cada nodo de elemento)
        xi = coordinates(elements(midpointload(i,5),1),1);
        yi = coordinates(elements(midpointload(i,5),1),2);
        xj = coordinates(elements(midpointload(i,5),2),1);
        yj = coordinates(elements(midpointload(i,5),2),2);
        xk = coordinates(elements(midpointload(i,5),3),1);
        yk = coordinates(elements(midpointload(i,5),3),2);
        xl = coordinates(elements(midpointload(i,5),4),1);
        yl = coordinates(elements(midpointload(i,5),4),2);
        
        l1= 1/2*sqrt((xj-xi)^2 + (yj-yi)^2);
        l2= 1/2*sqrt((xk-xj)^2 + (yk-yj)^2);
        A = l1*l2;
        
        %Familia de forma cuadrangular
        Ni = ((l1-x)*(l2-y))/(4*l1*l2);%symbolic
        Nj = ((l1+x)*(l2-y))/(4*l1*l2);
        Nk = ((l1+x)*(l2+y))/(4*l1*l2);
        Nl = ((l1-x)*(l2+y))/(4*l1*l2);
        
        N1= subs(Ni,[x,y],[midpointload(i,1),midpointload(i,2)]);
        N2= subs(Nj,[x,y],[midpointload(i,1),midpointload(i,2)]);
        N3= subs(Nk,[x,y],[midpointload(i,1),midpointload(i,2)]);
        N4= subs(Nl,[x,y],[midpointload(i,1),midpointload(i,2)]);
        
        %Para cada punto del elemento distribuyo la carga que se situa
        %entre medio del elemento
        ieqn = elements(midpointload(i,5),1)*2;         % Finds eq. number for the first node
        force(ieqn-1) = force(ieqn-1) + N1*midpointload(i,3)*A;   % add x force
        force(ieqn  ) = force(ieqn  ) + N1*midpointload(i,4)*A;   % add y force
        
        ieqn = elements(midpointload(i,5),2)*2;         % Finds eq. number for the second node
        force(ieqn-1) = force(ieqn-1) + N2*midpointload(i,3)*A;   % add x force
        force(ieqn  ) = force(ieqn  ) + N2*midpointload(i,4)*A;   % add y force
        
        ieqn = elements(midpointload(i,5),3)*2;         % Finds eq. number for the third node
        force(ieqn-1) = force(ieqn-1) + N3*midpointload(i,3)*A;   % add x force
        force(ieqn  ) = force(ieqn  ) + N3*midpointload(i,4)*A;   % add y force
        
        ieqn = elements(midpointload(i,5),4)*2;         % Finds eq. number for the fourth node
        force(ieqn-1) = force(ieqn-1) + N4*midpointload(i,3)*A;   % add x force
        force(ieqn  ) = force(ieqn  ) + N4*midpointload(i,4)*A;   % add y force
    end
end

%  Add point loads conditions to the force vector
for i = 1 : size(pointload,1)
    ieqn = (pointload(i,1)-1)*2 + pointload(i,2);       % Finds eq. number
    force(ieqn) = force(ieqn) + pointload(i,3);         % add the force
end

ttim = timing('Time for apply side and point load',ttim); %Reporting time

%  Applies the Dirichlet conditions and adjust the right hand side.
u = sparse ( nndof, 1 );
for i = 1 : size(fixnodes,1)
    ieqn = (fixnodes(i,1)-1)*2 + fixnodes(i,2);  %Finds the equation number
    u(ieqn) = fixnodes(i,3);                   %and store the solution in u
    fix(i) = ieqn;                         % and mark the eq as a fix value
end
force = force - StifMat * u;       % adjust the rhs with the known values

%  Compute the solution by solving StifMat * u = force for the
%  remaining unknown values of u.
FreeNodes = setdiff ( 1:nndof, fix ); % Finds the free node list and
% solve for it.
u(FreeNodes) = StifMat(FreeNodes,FreeNodes) \ force(FreeNodes);

%  Compute the reactions on the fixed nodes as a R = StifMat * u - F
reaction = sparse(nndof,1);
reaction(fix) = StifMat(fix,1:nndof) * u(1:nndof) - force(fix);

ttim = timing('Time  to solve the stifness matrix',ttim); %Reporting time

% Compute the stresses
Strnod = Stress(dmat,poiss,thick,pstrs,u);

ttim = timing('Time to  solve the  nodal stresses',ttim); %Reporting time

% Graphic representation.
ToGiD (['Problems/Results/', file_name],u,reaction,Strnod);

ttim = timing('Time  used to write  the  solution',ttim); %Reporting time
itim = toc;                                               %Close last tic
fprintf(1,'\n Total running time %12.6f \n',ttim);  %Reporting final time
