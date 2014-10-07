function [M] = TrStifCalmix ( h,node_i,node_j,nu_node )

%% TrStifCal Evaluates the mixload part of the stiffness matrix for a triangular element.
%
%  Parameters:
%
%    Input, nodes: contains the 2D coordinates nodes coordinates
%                  of the vertices.
%           h : pelicular coefficient

%    Output, M the mixload condition element matrix

l = sqrt((node_i(1)-node_j(1))^2 + (node_i(2)-node_j(2))^2);

%matrix that represents the mixsiload for this element using only form
%functions involved
switch nu_node
    case 1
        M = h*[0 0 0;0 l/3 l/6;0 l/6 l/3];
    case 2
        M = h*[l/3 0 l/6;0 0 0;l/6 0 l/3];
    case 3
        M = h*[l/3 l/6 0;l/6 l/3 0;0 0 0];
end
