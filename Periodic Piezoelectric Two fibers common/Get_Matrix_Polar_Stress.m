function [Smss,Smrr,Smrs]=Get_Matrix_Polar_Stress(R1,R2,x1,x2,a,Xce,Xdf,z)
[Sm22,Sm11,Sm12]=Get_Matrix_Stress(R1,R2,x1,x2,a,Xce,Xdf,z);

Smrr=(Sm22+Sm11)/2+(Sm11-Sm22)/2*cos(2*angle(z))+Sm12*sin(2*angle(z));
Smss=(Sm22+Sm11)/2-(Sm11-Sm22)/2*cos(2*angle(z))-Sm12*sin(2*angle(z));
Smrs=(Sm11-Sm22)/2*sin(2*angle(z))-Sm12*cos(2*angle(z));
end