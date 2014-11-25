Gid Post Results File 1.0 
### 
# MAT_FEM  V.1.0 
# 
Result "Displacements" "Load Analysis"  1  Vector OnNodes 
ComponentNames "X-Displ", "Y-Displ", "Z-Displ" 
Values 
     1   5.43083e-05         00000 
     2         00000         00000 
     3   6.41696e-04   1.42917e-06 
     4   6.65787e-04   1.14742e-04 
End Values 
# 
Result "Reaction" "Load Analysis"  1  Vector OnNodes 
ComponentNames "Rx", "Ry", "Rz" 
Values 
     1        00000 -2.45000e+05 
     2 -1.96000e+05  2.45000e+05 
     3        00000        00000 
     4        00000        00000 
End Values 
# 
Result "Stresses" "Load Analysis"  1  Matrix OnNodes 
ComponentNames "Sx", "Sy", "Sz", "Sxy", "Syz", "Sxz" 
Values 
     1 -1.21636e+05  2.15946e+04  1.26819e+05 
     2 -1.48728e+05 -4.31891e+04  2.36362e+05 
     3 -1.21636e+05  2.15946e+04  1.26819e+05 
     4 -9.45449e+04  8.63782e+04  1.72756e+04 
End Values 
 #