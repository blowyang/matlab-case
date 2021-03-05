function [delta1_Divided_a,delta2_Divided_a,delta_cof]=Get_Delta(R1,R2,x1,x2,a,Xcd1,Xcd2,Xe1,Xe2,S013_Divided_G0,S023_Divided_G0)
%delta_cof is equal module /G0
N=size(Xcd1,1)/4;
M=size(Xe1,1)/2;

S1=zeros(N,2*N);
S0=zeros(M,2*M);

for m=1:M
    S0(m,2*m-1)=1;
    S0(m,2*m)=1i; 
end
for n=1:N
    S1(n,2*n-1)=1;
    S1(n,2*n)=1i; 
end

% NP=100;
% del=a/(NP+1);

delta_cof=zeros(2,2);
%
% inte1=0;
% inte2=0;
% for np=1:NP
%     zk_AB=del*np-a/2+a*1i/2;    
%     temp1=Get_f0_Diff(a,R1,R2,x1,x2,Xcd1,Xe1,zk_AB);
%     temp2=Get_f0_Diff(a,R1,R2,x1,x2,Xcd2,Xe2,zk_AB);
%     inte1=inte1+real(temp1)*del;
%     inte2=inte2+real(temp2)*del;
% end
% delta_cof(2,1)=inte1;
% delta_cof(2,2)=inte2;
%
NP=100;
h=a/NP;
z0=-a/2+a/2*1i;
zN=a/2+a/2*1i;
temp1=0;
temp2=0;
temp1=temp1+h/6*Get_f0_Diff(a,R1,R2,x1,x2,Xcd1,Xe1,z0);
temp1=temp1+h/6*Get_f0_Diff(a,R1,R2,x1,x2,Xcd1,Xe1,zN);
temp2=temp2+h/6*Get_f0_Diff(a,R1,R2,x1,x2,Xcd2,Xe2,z0);
temp2=temp2+h/6*Get_f0_Diff(a,R1,R2,x1,x2,Xcd2,Xe2,zN);


for np=1:(NP-1)
    zn=-a/2+np*h+a/2*1i;
    temp1=temp1+h/3*Get_f0_Diff(a,R1,R2,x1,x2,Xcd1,Xe1,zn);
    temp2=temp2+h/3*Get_f0_Diff(a,R1,R2,x1,x2,Xcd2,Xe2,zn);
end

for np=1:NP
    zn=-a/2+(np-0.5)*h+a/2*1i;
    temp1=temp1+2*h/3*Get_f0_Diff(a,R1,R2,x1,x2,Xcd1,Xe1,zn);
    temp2=temp2+2*h/3*Get_f0_Diff(a,R1,R2,x1,x2,Xcd2,Xe2,zn);
end
delta_cof(2,1)=real(temp1);
delta_cof(2,2)=real(temp2);
%
NP=100;
h=a/NP;
z0=a/2-a/2*1i;
zN=a/2+a/2*1i;
temp1=0;
temp2=0;
temp1=temp1+h/6*Get_f0_Diff(a,R1,R2,x1,x2,Xcd1,Xe1,z0);
temp1=temp1+h/6*Get_f0_Diff(a,R1,R2,x1,x2,Xcd1,Xe1,zN);
temp2=temp2+h/6*Get_f0_Diff(a,R1,R2,x1,x2,Xcd2,Xe2,z0);
temp2=temp2+h/6*Get_f0_Diff(a,R1,R2,x1,x2,Xcd2,Xe2,zN);


for np=1:(NP-1)
    zn=a/2+(-a/2+np*h)*1i;
    temp1=temp1+h/3*Get_f0_Diff(a,R1,R2,x1,x2,Xcd1,Xe1,zn);
    temp2=temp2+h/3*Get_f0_Diff(a,R1,R2,x1,x2,Xcd2,Xe2,zn);
end

for np=1:NP
    zn=a/2+(-a/2+(np-0.5)*h)*1i;
    temp1=temp1+2*h/3*Get_f0_Diff(a,R1,R2,x1,x2,Xcd1,Xe1,zn);
    temp2=temp2+2*h/3*Get_f0_Diff(a,R1,R2,x1,x2,Xcd2,Xe2,zn);
end
delta_cof(1,1)=imag(temp1);
delta_cof(1,2)=imag(temp2);
%
% inte1=0;
% inte2=0;
% for np=1:NP
%     zk_DB=del*np*1i-a*1i/2+a/2;
%     temp1=Get_f0_Diff(a,R1,R2,x1,x2,Xcd1,Xe1,zk_DB);
%     temp2=Get_f0_Diff(a,R1,R2,x1,x2,Xcd2,Xe2,zk_DB);
%     inte1=inte1+imag(temp1)*del;
%     inte2=inte2+imag(temp2)*del;
% end
% delta_cof(1,1)=inte1;
% delta_cof(1,2)=inte2;
%

Stress_res(2,1)=S023_Divided_G0;
Stress_res(1,1)=S013_Divided_G0;
res=delta_cof\Stress_res;
delta1_Divided_a=res(1,1);
delta2_Divided_a=res(2,1);

end