function MisesMatrix=Get_MisesMatrix_Divided_G0(Xd,Xe,a,R1,Z)
%��û����ϵĻ��������Ӧ��
S013_Stress_Divided_G0=Get_S013_Stress_Divided_G0(Xd,Xe,a,R1,Z);
S023_Stress_Divided_G0=Get_S023_Stress_Divided_G0(Xd,Xe,a,R1,Z);

MisesMatrix=(3*(S013_Stress_Divided_G0^2+S023_Stress_Divided_G0^2))^0.5;
end