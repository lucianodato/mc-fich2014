%=======================================================================
% MAT-fem 1.0  - MAT-fem is a learning tool for undestanding 
%                the Finite Element Method with MATLAB and GiD
%=======================================================================
% PROBLEM TITLE = Titulo del problema
%
%  Material Properties
%
  young =   2.1e9;
  poiss =              0.30000 ;
  denss = 2400.00000 ;
  pstrs =  1 ;
  thick =              0.00500 ;
%
% Coordinates
%
global coordinates
coordinates = [
         0.00000   ,         0.00000  ;
         1.00000   ,         0.00000  ;
         2.00000   ,         0.00000  ;
         0.00000   ,         1.00000  ;
         0.50000   ,         1.00000  ;
         1.50000   ,         1.00000  ;
         2.00000   ,         1.00000  ;
        ] ; 
%
% Elements
%
global elements
elements = [
      1   ,      5   ,      4   ; 
      1   ,      2   ,      5   ; 
      2   ,      6   ,      5   ; 
      2   ,      3   ,      6   ; 
      3   ,      7   ,      6  ] ; 
%
% Fixed Nodes
%
fixnodes = [
      1  , 1 ,    0.00000  ;
      1  , 2 ,    0.00000  ;
      4  , 1 ,    0.00000  ;
      4  , 2 ,    0.00000  ;
      3  , 2 ,    0.00000  ] ;
%
% Point loads
%
pointload = [
     7  , 1 ,    3000.00000  ;
     7  , 2 ,   5000.00000  ] ;
%
%Uniform Side loads
%
sideload = [6  , 5 , 0.00000 ,-2000.00000];

