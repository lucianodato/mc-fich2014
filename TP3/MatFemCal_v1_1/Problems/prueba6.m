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
         2.50000   ,         5.00000  ;
         0.00000   ,         2.50000  ;
         2.50000   ,         2.50000  ;
         0.00000   ,         0.00000  ;
         5.00000   ,         5.00000  ;
         2.50000   ,         0.00000  ;
         5.00000   ,         2.50000  ;
         5.00000   ,         0.00000  ] ; 
%
% Elements
%
global elements
elements = [
      7   ,      4   ,      5   ; 
      4   ,      3   ,      5   ; 
      9   ,      8   ,      7   ; 
      8   ,      4   ,      7   ; 
      4   ,      2   ,      3   ; 
      2   ,      1   ,      3   ; 
      8   ,      6   ,      4   ; 
      6   ,      2   ,      4   ] ; 
%
% Fixed Nodes
%
fixnodes = [
      5  ,     0.00000  ;
      6  ,   100.00000  ;
      7  ,     0.00000  ;
      8  ,   100.00000  ;
      9  ,   100.00000  ] ;
%
% Punctual Fluxes
%
pointload = [ ] ;
%
% Side loads
%
sideload = [
      3  ,      5  ,    2.00000   ;
      1  ,      3  ,    2.00000  ];

