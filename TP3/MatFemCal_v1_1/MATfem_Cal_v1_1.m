%% MAT-femCal
%
% Clear memory and variables.
clear

% The variables are readed as a MAT-femCal subroutine
% kx   = Heat transfer coeffcient in x direction
% kx   = Heat transfer coeffcient in x direction
% heat = Heat source per area unit
% coordinates = [ x , y ] coordinate matrix nnode x ndime (2)
% elements    = [ inode, jnode, knode ] element conectivities matrix
%               nelem x nnode; nnode = 3 for triangular elements and
%               nnode = 4 for cuadrilateral elements
% fixnodes    = [node number, fixed value] matrix with
%               restricted temperatures restrictions.
% pointload   = [node number load value] matrix with
%               flux load.
% SideLoad    = [node number i, node number j, normal flux] matrix with
%               line definition and uniform normal flux applied
% MixLoad    = [node number i, node number j, normal flux,fixed value] 
%               matrix with line definition and uniform normal flux 
%               applied and restriced temperature restrictions
% MixLoadelements = elements that have mixload conditions (same order as mixload)
%Mixed Load Conditions initialization %@ Agregado
h=1;
mixloadelements = [];
mixload = [];%Manual completition

%Transient Flag %@ Agregado
transient = 0;

if (transient == 1) %@ Agregado
    % VARIABLES DE ENTRADA INICIAL
    
    rho = 7897;
    cp = 0.108;
    tmax = 100000;
    dt = 2000;
    contador = 0;
    paso_graph = floor((tmax/dt+1)/10);
    
    flag_euler = 1; % 0 forward, 1 backward, 2 crank-nicholson
end 
% FIN VARIABLES
%@ Agregado

file_name = input('Enter the file name :','s');

tic;                   % Start clock
ttim = 0;              % Initialize time counter
run(eval (file_name));      % Read input file

% Finds basics dimentions
npnod  = size(coordinates,1);      % Number of nodes
nndof  = npnod;                    % Number of total DOF
nelem  = size(elements,1);         % Number of elements
nnode  = size(elements,2);         % Number of nodes per element
neleq  = nnode;                    % Number of DOF per element

ttim = timingcal('Time needed to read the input file',ttim); %Reporting time

% Dimension the global matrices.
StifMat = sparse ( nndof , nndof );  % Create the global stiffness matrix
if (transient == 1)
    C_Mat = sparse ( nndof , nndof );    %@ AGREGADO
end
force   = sparse ( nndof , 1 );      % Create the global force vector

%  Material properties (Constant over the domain).
dmat = consttCal(kx,ky);

ttim = timingcal('Time needed to  set initial values',ttim); %Reporting time

%  Element cycle.
for ielem = 1 : nelem
    
    % Recover element properties
    lnods = elements(ielem,:);                        % Elem. conectivities
    coord(1:nnode,:) = coordinates(lnods(1:nnode),:); % Elem. coordinates
    
    % Evaluates the elemental stiffnes matrix and mass force vector.
    if (nnode == 3)
        
        [ElemMat,ElemFor] = TrStifCal(coord,dmat,heat); % 3 Nds Triangle
        
        if (find(mixloadelements == ielem))
            index = find(mixloadelements == ielem);%index of the element in mixload
            node_i= coordinates(mixload(index,1));%border node of the element
            node_j= coordinates(mixload(index,2));%border node
            
            %We should know also the node that is not used in order to
            %build the elemental mixload matrix
            if(lnods(1) ~= mixload(index,1) && lnods(1) ~= mixload(index,2))
                nu_node = 1;
            end
            if(lnods(2) ~= mixload(index,1) && lnods(2) ~= mixload(index,2))
                nu_node = 2;
            end
            if(lnods(3) ~= mixload(index,1) && lnods(3) ~= mixload(index,2))
                nu_node = 3;
            end
            
            ElemMat = ElemMat + TrStifCalmix(h,node_i,node_j,nu_node); %mixload part of the stiffness matrix is added
        end
        
        if (transient == 1)
            [Elem_C_Mat] = obtener_C_mat_tri(coord, rho, cp);% @Agregado
        end
    else
        [ElemMat,ElemFor] = QdStifCal(coord,dmat,heat); % 4 Nds Quad.
        
        %falta condiciones mixtas para cuadrangulos
        
        if (transient == 1)
            [Elem_C_Mat] = obtener_C_mat_quad(coord, rho, cp);% @Agregado
        end
    end
    
    % Finds the equation number list for the i-th element
    eqnum = [];                                  % Clear the list
    for i =1 : nnode                             % Node cicle
        eqnum = [eqnum,lnods];                     % Build the equation
    end                                          % number list
    
    % Assamble the force vector and the stiffnes matrix
    for i = 1 : neleq
        force(eqnum(i)) = force(eqnum(i)) + ElemFor(i);
        for j = 1 : neleq
            
            StifMat(eqnum(i),eqnum(j)) = StifMat(eqnum(i),eqnum(j)) + ...
                ElemMat(i,j);
            if (transient == 1)
                C_Mat(eqnum(i),eqnum(j)) = C_Mat(eqnum(i),eqnum(j)) + ...     %@ Agregado
                    Elem_C_Mat(i,j);                                          %@ Agregado
            end
        end
    end
    
end  % End element cicle

ttim = timingcal('Time to assamble the global system',ttim); %Reporting time

%  Add mix side forces to the force vector %@ Agregado
for i = 1 : size(mixload,1)
    x=coordinates(mixload(i,1),:)-coordinates(mixload(i,2),:);
    l = sqrt(x*transpose(x));       % Finds the lenght of the side
    ieqn = mixload(i,1);           % Finds eq. number for the first node
    force(ieqn) = force(ieqn) + l*mixload(i,3)/2 + h*(l/2)*mixload(i,4);% add puntual heat
    
    ieqn = mixload(i,2);            % Finds eq. number for the second node
    force(ieqn) = force(ieqn) + l*mixload(i,3)/2 + h*(l/2)*mixload(i,4);% add puntual heat
end %@ Agregado

%  Add side forces to the force vector
for i = 1 : size(sideload,1)
    x=coordinates(sideload(i,1),:)-coordinates(sideload(i,2),:);
    l = sqrt(x*transpose(x));       % Finds the lenght of the side
    ieqn = sideload(i,1);           % Finds eq. number for the first node
    force(ieqn) = force(ieqn) + l*sideload(i,3)/2;     % add puntual heat
    
    ieqn = sideload(i,2);            % Finds eq. number for the second node
    force(ieqn) = force(ieqn) + l*sideload(i,3)/2;      % add puntual heat
end

%  Add point loads conditions to the force vector
for i = 1 : size(pointload,1)
    ieqn = pointload(i,1);                              % Finds eq. number
    force(ieqn) = force(ieqn) + pointload(i,2);         % add the force
end

ttim = timingcal('Time for apply side and point load',ttim); %Reporting time

%  Applies the Dirichlet conditions and adjust the right hand side.
u = sparse ( nndof, 1 );
for i = 1 : size(fixnodes,1)
    ieqn = fixnodes(i,1);                      %Finds the equation number
    u(ieqn) = fixnodes(i,2);                   %and store the solution in u
    fix(i) = ieqn;                         % and mark the eq as a fix value
end
force = force - StifMat * u;       % adjust the rhs with the known values



%  Compute the solution by solving StifMat * u = force for the
%  remaining unknown values of u.

FreeNodes = setdiff ( 1:nndof, fix ); % Finds the free node list and
% solve for it.
%If it's a transient problem
if (transient == 1)
    for t = 0 : dt : tmax
        
        if (t == 0)
            u_anterior = zeros(size(force));
        else
            u_anterior = u;
        end
        
        %@ AGREGADO
        u(FreeNodes) = euler(StifMat(FreeNodes, FreeNodes), ...
            C_Mat(FreeNodes, FreeNodes), ...
            force(FreeNodes), ...
            dt, ...
            u_anterior(FreeNodes), ...
            flag_euler);
        
        %  Compute the reactions on the fixed nodes as a R = StifMat * u - F
        reaction = sparse(nndof,1);
        
        if (flag_euler) % backward
            reaction(fix) = StifMat(fix,1:nndof) * u(1:nndof) + ...
                + C_Mat(fix, 1:nndof) * (u(1:nndof) - u_anterior(1:nndof)) - force(fix);
            
        else % forward
            reaction(fix) = StifMat(fix,1:nndof) * u_anterior(1:nndof) + ...
                + C_Mat(fix, 1:nndof) * (u(1:nndof) - u_anterior(1:nndof)) - force(fix);
        end
        ttim = timingcal('Time  to solve the stifness matrix',ttim); %Reporting time
        
        % Compute the stresses
        Strnod = StressCal(dmat,u);
        
        if (mod(t/dt, paso_graph) == 0)
            contador = contador + 1;
            
            ToGiDCal ([ 'results/', file_name, '.', num2str(contador) ],u,reaction,Strnod);
        end
        
    end
    
else
    %if isn't a transient problem
    u(FreeNodes) = StifMat(FreeNodes,FreeNodes) \ force(FreeNodes);
    
    %  Compute the reactions on the fixed nodes as a R = StifMat * u - F
    reaction = sparse(nndof,1);
    reaction(fix) = StifMat(fix,1:nndof) * u(1:nndof) - force(fix);
    
    ttim = timingCal('Time  to solve the stifness matrix',ttim); %Reporting time
    
    % Compute the stresses
    Strnod = StressCal(dmat,u);
    
    ttim = timingCal('Time to  solve the  nodal stresses',ttim); %Reporting time
    
    % Graphic representation.
    ToGiDCal (file_name,u,reaction,Strnod);
end

ttim = timingcal('Time  used to write  the  solution',ttim); %Reporting time
itim = toc;                                               %Close last tic
fprintf(1,'\n Total running time %12.6f \n',ttim);  %Reporting final time
