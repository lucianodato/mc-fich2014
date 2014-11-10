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
 heat=              100.00000 ;
%
% Coordinates
%
global coordinates
coordinates = [
         0.00000   ,         0.00000  ;
         1.00000   ,         0.00000  ;
         1.00000   ,         1.00000  ;
         0.50000   ,         1.30000  ;
         0.00000   ,         1.00000  ] ; 
%
% Elements
%
global elements
elements = [
      1   ,      2   ,      5   ; 
      2  ,      3   ,      5   ;
      3   ,      4   ,      5   ] ; 
%
% Fixed Nodes
%
fixnodes = [
      1  ,    100.00000 ;
      2  ,    100.00000 ] ;
%
% Punctual Fluxes
%
pointload = [ ] ;
%
% Side loads
%
sideload = [
      1  ,      5  ,    0.00000   ;
      3  ,      4  ,    0.00000   ;
      4  ,      5  ,    20.00000  ;
      2  ,      3  ,   -20.00000   ];

