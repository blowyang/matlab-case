function MisesMatrix=Get_MisesMatrix(a,R,Xce,Xdf,z)
[Sm22,Sm11,Sm12]=Get_Matrix_Stress(a,R,Xce,Xdf,z);
%��û����ϵĻ��������Ӧ��

MisesMatrix=((Sm22^2+Sm11^2+(Sm22-Sm11)^2+6*Sm12^2)/2)^0.5;

end