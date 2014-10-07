Gid Post Results File 1.0 
### 
# MAT_FEM  V.1.0 
# 
Result "Temperatures" "Load Analysis"  1  Scalar OnNodes 
ComponentNames "Temperature" 
Values 
     1         00050 
     2   2.60000e+01 
     3         00000 
     4         00000 
End Values 
# 
Result "ReactiveFluxes" "Load Analysis"  1  Scalar OnNodes 
ComponentNames "Fx" 
Values 
     1  8.60000e+01 
     2        00000 
     3 -5.00000e+01 
     4 -1.30000e+01 
End Values 
# 
Result "Fluxes" "Load Analysis"  1  Vector OnNodes 
ComponentNames "Fx", "Fy", "Fz" 
Values 
     1       -00050  2.40000e+01  0.0 
     2 -3.80000e+01  1.20000e+01  0.0 
     3 -3.80000e+01  1.20000e+01  0.0 
     4 -2.60000e+01        00000  0.0 
End Values 
