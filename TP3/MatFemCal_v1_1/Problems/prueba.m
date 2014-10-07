%=======================================================================
% MAT-femcCal 1.0  - MAT-femCal is a learning tool for undestanding 
%                    the Finite Element Method with MATLAB and GiD
%=======================================================================
% PROBLEM TITLE = Titulo del problema
%
%  Material Properties
%
  kx =              1.00000 ;
  ky =              1.00000 ;
 heat=              0.00000 ;
%
% Coordinates
%
global coordinates
coordinates = [
         0.00000   ,         1.00000  ;
         0.00000   ,         0.00000  ;
         1.00000   ,         1.00000  ;
         1.00000   ,         0.00000  ] ; 
%
% Elements
%
global elements
elements = [
      2   ,      4   ,      3   ; 
      3   ,      1   ,      2   ] ; 
%
% Fixed Nodes
%
fixnodes = [
      1  ,    50.00000  ;
      3  ,    0.00000  ; 
      4  ,    0.00000] ;
%
% Punctual Fluxes
%
pointload = [ ] ;
%
% Side loads
%
sideload = [
      2  ,      4  ,    0.00000   ;
      1  ,      2  ,    2.00000  ];

