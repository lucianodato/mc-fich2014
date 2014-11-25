Gid Post Results File 1.0 
### 
# MAT_FEM  V.1.0 
# 
Result "Displacements" "Load Analysis"  1  Vector OnNodes 
ComponentNames "X-Displ", "Y-Displ", "Z-Displ" 
Values 
     1         00000         00000 
     2   3.21617e-04   1.50775e-04 
     3   4.82155e-04         00000 
     4         00000         00000 
     5   1.12286e-04  -7.15539e-05 
     6   4.38125e-04   2.04558e-04 
     7   4.62285e-04   8.41159e-04 
End Values 
# 
Result "Reaction" "Load Analysis"  1  Vector OnNodes 
ComponentNames "Rx", "Ry", "Rz" 
Values 
     1 -1.55991e+03 -2.52264e+02 
     2        00000        00000 
     3        00000 -3.20804e+03 
     4 -1.44009e+03  4.84308e+02 
     5        00000        00000 
     6        00000        00000 
     7        00000        00000 
End Values 
# 
Result "Stresses" "Load Analysis"  1  Matrix OnNodes 
ComponentNames "Sx", "Sy", "Sz", "Sxy", "Syz", "Sxz" 
Values 
     1  5.79353e+05  1.95172e+04 -1.64993e+04 
     2  5.97985e+05  1.07964e+05  5.85351e+04 
     3  5.76872e+05  1.19191e+06  4.59898e+05 
     4  5.18242e+05  1.55472e+05 -1.15587e+05 
     5  6.17433e+05  2.33796e+04  5.08428e+04 
     6  6.15779e+05  8.04973e+05  3.68441e+05 
     7  6.93847e+05  1.97459e+06  1.01231e+06 
End Values 
 #