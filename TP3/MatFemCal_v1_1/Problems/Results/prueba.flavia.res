Gid Post Results File 1.0 
### 
# MAT_FEM  V.1.0 
# 
Result "Temperatures" "Load Analysis"  1  Scalar OnNodes 
ComponentNames "Temperature" 
Values 
     1         00050 
     2         00052 
     3         00050 
     4         00052 
End Values 
# 
Result "ReactiveFluxes" "Load Analysis"  1  Scalar OnNodes 
ComponentNames "Fx" 
Values 
     1        00102 
     2        00000 
     3        00102 
     4        00000 
End Values 
# 
Result "Fluxes" "Load Analysis"  1  Vector OnNodes 
ComponentNames "Fx", "Fy", "Fz" 
Values 
     1        00000       -00002  0.0 
     2        00000       -00002  0.0 
     3        00000       -00002  0.0 
     4        00000       -00002  0.0 
End Values 
