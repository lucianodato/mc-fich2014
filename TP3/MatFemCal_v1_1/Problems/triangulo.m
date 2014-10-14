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
         5.00000   ,         5.00000  ;
         0.00000   ,         0.00000  ;
         5.00000   ,         0.00000  ] ; 
%
% Elements
%
global elements
elements = [
      3   ,      4   ,      1   ; 
      4   ,      2   ,      1   ] ; 
%
% Fixed Nodes
%
fixnodes = [
	  2  ,   100.00000  ;
      3  ,   0.00000  ;
      4  ,   0.00000  ] ;
%
% Punctual Fluxes
%
pointload = [ ] ;
%
% Side loads
%sideload = [];
sideload = [
      1  ,      3  ,    2.00000  ];

