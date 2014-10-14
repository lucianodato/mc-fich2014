%=======================================================================
% MAT-fem 1.0  - MAT-fem is a learning tool for undestanding 
%                the Finite Element Method with MATLAB and GiD
%=======================================================================
% PROBLEM TITLE = Titulo del problema
%
%  Material Properties
%
  young =   200000000000.00000 ;
  poiss =              0.30000 ;
  denss = 0.00 ;
  pstrs =  1 ;
  thick =              0.10000 ;
%
% Coordinates
%
global coordinates
coordinates = [
         0.00000   ,         0.28000  ;
         0.00000   ,         0.14000  ;
         0.20000   ,         0.28000  ;
         0.20000   ,         0.14000  ;
         0.00000   ,         0.00000  ;
         0.20000   ,         0.00000  ;
         0.40000   ,         0.28000  ;
         0.40000   ,         0.14000  ;
         0.40000   ,         0.00000  ;
         0.60000   ,         0.28000  ;
         0.60000   ,         0.14000  ;
         0.60000   ,         0.00000  ;
         0.80000   ,         0.28000  ;
         0.80000   ,         0.14000  ;
         0.80000   ,         0.00000  ] ; 
%
% Elements
%
global elements
elements = [
      6   ,      4   ,      2   ,      5   ; 
      9   ,      8   ,      4   ,      6   ; 
     12   ,     11   ,      8   ,      9   ; 
     15   ,     14   ,     11   ,     12   ; 
      4   ,      3   ,      1   ,      2   ; 
      8   ,      7   ,      3   ,      4   ; 
     11   ,     10   ,      7   ,      8   ; 
     14   ,     13   ,     10   ,     11   ] ; 
%
% Fixed Nodes
%
fixnodes = [
      5  , 1 ,    0.00000  ;
      5  , 2 ,    0.00000  ;
     15  , 2 ,    0.00000  ] ;
%
% Point loads
%
pointload = [
     10  , 1 , 3500.00000  ;
     10  , 2 , 6062.17000  ] ;
%
%Uniform Side loads
%
sideload = [ ];

