clear;
clc;

R2=1;
VF=0.15;
a=(pi/VF)^0.5*R2;
G0=1250;
gama1=0.1;
G1=G0/gama1;
gama2=0.1;
G2=G1/gama2;
R1=R2+0.001;
lambda=R2/R1;
lambda1=R1/a ;

N=30;
M=25;
K=20;
[E_Matrix,F_Matrix,G_Matrix,N_Matrix]=Get_E_N_Matrix(N,lambda,gama1,gama2);
B_Matrix=Get_B_Matrix(N,M,gama1,lambda1);

C_Matrix=Get_C_Matrix(K,N,a,R1);
D_Matrix=Get_D_Matrix(K,M,a);

H_Matrix=Get_H_Matrix(D_Matrix,C_Matrix,N_Matrix,B_Matrix);

% P1=zeros(4*K,1);
% P2=zeros(4*K,1);
% 
% for i=1:K
%     P1(i,1)=1;
%     P2(i+K,1)=1;
% end
% 
% %Çódelt
% Xe1=(H_Matrix'*H_Matrix)\(H_Matrix')*P1;
% Xe2=(H_Matrix'*H_Matrix)\(H_Matrix')*P2;
% 
% Xd1=-N_Matrix*B_Matrix*Xe1;
% Xd2=-N_Matrix*B_Matrix*Xe2;
[Xa1,Xa2,Xd1,Xd2,Xe1,Xe2]=Get_Xa_Xe_SubElement(K,H_Matrix,B_Matrix,E_Matrix,N_Matrix);
% NP=100;
% del=a/(NP+1);
% 
% S1=zeros(N,2*N);
% S2=zeros(M,2*M);
% for n=1:N
%     S1(n,2*n-1)=1;
%     S1(n,2*n)=1i; 
% end
% for m=1:M
%     S2(m,2*m-1)=1;
%     S2(m,2*m)=1i; 
% end
% integral1=0;
% integral2=0;
% x_cof=zeros(2,2);
% for np=1:NP
%     z_AB=del*np-a/2+a*1i/2;
%     L1=zeros(N,1);
%     L2=zeros(M,1);
%     for n=1:N
%         L1(n,1)=n*(R1/z_AB)^(n+1);
%     end
%     for m=1:M
%         L2(m,1)=Get_P_Diff(m,a,z_AB);
%     end    
%     integral1=integral1+real(Xe2.'*S2.'*L2-1/R1*Xd2.'*S1.'*L1)*del;
%     integral2=integral2+real(Xe1.'*S2.'*L2-1/R1*Xd1.'*S1.'*L1)*del;
% end
% x_cof(1,1)=integral1;
% x_cof(1,2)=integral2;
% integral1=0;
% integral2=0;
% for np=1:NP
%     z_DB=del*np*1i-a*1i/2+a/2;
%     L1=zeros(N,1);
%     L2=zeros(M,1);
%     for n=1:N
%         L1(n,1)=n*(R1/z_DB)^(n+1);
%     end
%     for m=1:M
%         L2(m,1)=Get_P_Diff(m,a,z_DB);
%     end    
%     integral1=integral1+imag(Xe2.'*S2.'*L2-1/R1*Xd2.'*S1.'*L1)*del;
%     integral2=integral2+imag(Xe1.'*S2.'*L2-1/R1*Xd1.'*S1.'*L1)*del;
% end
% x_cof(2,1)=integral1;
% x_cof(2,2)=integral2;
% y_res(1,1)=0.1;
% y_res(2,1)=0.1;
% res=x_cof\y_res;
S013_Divided_G0=0.1;
S023_Divided_G0=0.1;
[delta1_Divided_a,delta2_Divided_a,delta_cof]=Get_Delta(R1,a,Xd1,Xd2,Xe1,Xe2,S013_Divided_G0,S023_Divided_G0);
delta1=delta1_Divided_a*a;
delta2=delta2_Divided_a*a;


% delta1=1;
% delta2=0;
% Xd=delta1*Xd2+delta2*Xd1;
% Xe=delta1*Xe2+delta2*Xe1;
% NP=100;
% del=a/(NP+1);
% w_AB=zeros(NP,1);
% w_CD=zeros(NP,1);
% w_AB_Sub_CD=zeros(NP,1);
% w_DB=zeros(NP,1);
% w_CA=zeros(NP,1);
% w_DB_Sub_CA=zeros(NP,1);
% S23_AB=zeros(NP,1);
% S13_DB=zeros(NP,1);
% S13_CA=zeros(NP,1);
% S23_CD=zeros(NP,1);
% S13_DB_Sub_CA=zeros(NP,1);
% S23_AB_Sub_CD=zeros(NP,1);
% for np=1:NP
%     z_AB=del*np-a/2+a*1i/2;
%     z_DB=del*np*1i-a*1i/2+a/2;
%     z_CD=conj(z_AB);
%     z_CA=conj(-z_DB);
%     
%     w_AB(np,1)=Get_W0_Dispalcement(Xd,Xe,a,R1,z_AB); 
%     w_CD(np,1)=Get_W0_Dispalcement(Xd,Xe,a,R1,z_CD);
%     w_CA(np,1)=Get_W0_Dispalcement(Xd,Xe,a,R1,z_CA); 
%     w_DB(np,1)=Get_W0_Dispalcement(Xd,Xe,a,R1,z_DB);
%     
%     w_AB_Sub_CD(np,1)=(w_AB(np,1)-w_CD(np,1))/(delta2);
%     w_DB_Sub_CA(np,1)=(w_DB(np,1)-w_CA(np,1))/(delta1);
%     
%     S23_AB(np,1)=Get_S023_Stress_Divided_G0(Xd,Xe,a,R1,z_AB);
%     S23_CD(np,1)=Get_S023_Stress_Divided_G0(Xd,Xe,a,R1,z_CD);
%     S13_DB(np,1)=Get_S013_Stress_Divided_G0(Xd,Xe,a,R1,z_DB);
%     S13_CA(np,1)=Get_S013_Stress_Divided_G0(Xd,Xe,a,R1,z_CA);
%  
%     S13_DB_Sub_CA(np,1)=S13_DB(np,1)-S13_CA(np,1);
%     S23_AB_Sub_CD(np,1)=S23_AB(np,1)-S23_CD(np,1);
% end
