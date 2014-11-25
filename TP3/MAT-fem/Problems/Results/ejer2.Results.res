Gid Post Results File 1.0 
### 
# MAT_FEM  V.1.0 
# 
Result "Displacements" "Load Analysis"  1  Vector OnNodes 
ComponentNames "X-Displ", "Y-Displ", "Z-Displ" 
Values 
     1   4.57169e-05         00000 
     2         00000         00000 
     3   6.41543e-04  -6.91990e-05 
     4   6.61666e-04   4.66647e-05 
End Values 
# 
Result "Reaction" "Load Analysis"  1  Vector OnNodes 
ComponentNames "Rx", "Ry", "Rz" 
Values 
     1        00000 -1.87433e+06 
     2 -1.96000e+06  2.81633e+06 
     3        00000        00000 
     4        00000        00000 
End Values 
# 
Result "Stresses" "Load Analysis"  1  Matrix OnNodes 
ComponentNames "Sx", "Sy", "Sz", "Sxy", "Syz", "Sxz" 
Values 
     1 -1.21796e+05 -4.78058e+04  1.26022e+05 
     2 -1.48409e+05 -1.13722e+05  2.37956e+05 
     3 -1.21796e+05 -4.78058e+04  1.26022e+05 
     4 -9.51823e+04  1.81100e+04  1.40887e+04 
End Values 
 #