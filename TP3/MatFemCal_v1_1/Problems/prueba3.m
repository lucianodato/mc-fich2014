%=======================================================================
% MAT-femcCal 1.0  - MAT-femCal is a learning tool for undestanding 
%                    the Finite Element Method with MATLAB and GiD
%=======================================================================
% PROBLEM TITLE = Titulo del problema
%
%  Material Properties
%
  kx =              2.00000 ;
  ky =              2.00000 ;
 heat=              0.00000 ;
%
% Coordinates
%
global coordinates
coordinates = [
         0.00000   ,         5.00000  ;
         0.00000   ,         0.00000  ;
         5.00000   ,         5.00000  ;
         5.00000   ,         0.00000  ] ; 
%
% Elements
%
global elements
elements = [
      2   ,      1   ,      4   ; 
      3   ,      1   ,      4  ] ; 
%
% Fixed Nodes
%
fixnodes = [
      2  ,     0.00000  ;
      3  ,   100.00000  ;
      4  ,     0.00000  ] ;
%
% Punctual Fluxes
%
pointload = [ ] ;
%
% Side loads
%
sideload = [
      1  ,      2  ,    2.00000  ];

