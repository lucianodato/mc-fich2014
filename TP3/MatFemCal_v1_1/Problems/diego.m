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
         0.50000   ,         1.00000  ;
         0.00000   ,         0.50000  ;
         0.50000   ,         0.50000  ;
         0.00000   ,         0.00000  ;
         1.00000   ,         1.00000  ;
         0.50000   ,         0.00000  ;
         1.00000   ,         0.50000  ;
         1.00000   ,         0.00000  ] ; 
%
% Elements
%
global elements
elements = [
      7   ,      4   ,      3   ,      5   ; 
      9   ,      8   ,      4   ,      7   ; 
      4   ,      2   ,      1   ,      3   ; 
      8   ,      6   ,      2   ,      4   ] ; 
%
% Fixed Nodes
%
fixnodes = [
      1  ,    10.00000  ;
      3  ,    10.00000  ;
      5  ,    10.00000  ;
      6  ,   100.00000  ;
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
      5  ,      7  ,    0.00000   ;
      7  ,      9  ,    0.00000   ;
      2  ,      1  ,    0.00000   ;
      6  ,      2  ,    0.00000  ];

