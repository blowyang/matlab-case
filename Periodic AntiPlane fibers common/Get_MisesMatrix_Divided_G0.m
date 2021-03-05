function MisesMatrix=Get_MisesMatrix_Divided_G0(a,R1,R2,x1,x2,Xcd,Xe,z)
%获得基体上的环向和轴向应力
S013_Stress_Divided_G0=Get_S013_Stress_Divided_G0(a,R1,R2,x1,x2,Xcd,Xe,z);
S023_Stress_Divided_G0=Get_S023_Stress_Divided_G0(a,R1,R2,x1,x2,Xcd,Xe,z);

MisesMatrix=(3*(S013_Stress_Divided_G0^2+S023_Stress_Divided_G0^2))^0.5;
end