function [Sigma_zs_Divided_G0,Sigma_zr_Divided_G0]=Get_SigmaMatrix_Divided_G0(Xd,Xe,a,R1,Z)
%获得基体上的环向和轴向应力
S013_Stress_Divided_G0=Get_S013_Stress_Divided_G0(Xd,Xe,a,R1,Z);
S023_Stress_Divided_G0=Get_S023_Stress_Divided_G0(Xd,Xe,a,R1,Z);

Sigma_zs_Divided_G0=real(exp(1i*angle(Z))*(S023_Stress_Divided_G0+S013_Stress_Divided_G0*1i));
Sigma_zr_Divided_G0=imag(exp(1i*angle(Z))*(S023_Stress_Divided_G0+S013_Stress_Divided_G0*1i));
end