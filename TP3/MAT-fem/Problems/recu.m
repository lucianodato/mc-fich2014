%=======================================================================
% MAT-fem 1.0  - MAT-fem is a learning tool for undestanding 
%                the Finite Element Method with MATLAB and GiD
%=======================================================================
% PROBLEM TITLE = Titulo del problema
%
%  Material Properties
%
  young =   20e9 ;
  poiss =              0.30000 ;
  denss = 0.00 ;
  pstrs =  1 ;
  thick =              0.10000 ;
%
% Coordinates
%
global coordinates
coordinates = [
         0.00000   ,         0.00000  ;
         8.00000   ,         0.00000  ;
         4.00000   ,         20.00000  ;
         0.00000   ,         20.00000 
         ] ; 
%
% Elements
%
global elements
elements = [
      1   ,      3   ,      4   ; 
      1   ,      2   ,      3    
      ] ; 
%
% Fixed Nodes
%
fixnodes = [
      2  , 1 ,    0.00000  ;
      2  , 2 ,    0.00000  ;
      1	, 2 ,    0.00000  ] ;
%
% Point loads
%
pointload = [ ] ;
%
%Uniform Side loads
%
sideload = [1, 4,9800.00000,0.00000 ];

