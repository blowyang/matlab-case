function S013_Stress_Divided_G0=Get_S013_Stress_Divided_G0(a,R1,R2,x1,x2,Xcd,Xe,z)

f0_Diff=Get_f0_Diff(a,R1,R2,x1,x2,Xcd,Xe,z);

S013_Stress_Divided_G0=imag(f0_Diff);
% N=size(Xcd,1)/4;
% M=size(Xe,1)/2;
% S1=zeros(N,2*N);
% S0=zeros(M,2*M);
% for n=1:N
%     S1(n,2*n-1)=1;
%     S1(n,2*n)=1i; 
% end
% for m=1:M
%     S0(m,2*m-1)=1;
%     S0(m,2*m)=1i; 
% end
% 
% L4=zeros(N,1);
% L5=zeros(N,1);
% L6=zeros(M,1);
% 
% for n=1:N
%     L4(n,1)=-n/R1*(R1/(z-x1))^(n+1);
%     L5(n,1)=-n/R2*(R2/(z-x2))^(n+1);
% end
% for m=1:M
%     L6(m,1)=Get_P_Diff(m,a,z);
% end 
% S013_Stress_Divided_G0=imag([L4.'*S1,L5.'*S1]*Xcd+L6.'*S0*Xe);
end