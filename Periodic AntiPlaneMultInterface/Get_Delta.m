function [delta1_Divided_a,delta2_Divided_a,delta_cof]=Get_Delta(R1,a,Xd1,Xd2,Xe1,Xe2,S013_Divided_G0,S023_Divided_G0)
%delta_cof is equal module /G0
N=size(Xd1,1)/2;
M=size(Xe1,1)/2;

NP=100;
del=a/(NP+1);

S1=zeros(N,2*N);
S2=zeros(M,2*M);
for n=1:N
    S1(n,2*n-1)=1;
    S1(n,2*n)=1i; 
end
for m=1:M
    S2(m,2*m-1)=1;
    S2(m,2*m)=1i; 
end
integral1=0;
integral2=0;
delta_cof=zeros(2,2);
for np=1:NP
    z_AB=del*np-a/2+a*1i/2;
    L1=zeros(N,1);
    L2=zeros(M,1);
    for n=1:N
        L1(n,1)=n*(R1/z_AB)^(n+1);
    end
    for m=1:M
        L2(m,1)=Get_P_Diff(m,a,z_AB);
    end    
    integral1=integral1+real(Xe2.'*S2.'*L2-1/R1*Xd2.'*S1.'*L1)*del;
    integral2=integral2+real(Xe1.'*S2.'*L2-1/R1*Xd1.'*S1.'*L1)*del;
end
delta_cof(1,1)=integral1;
delta_cof(1,2)=integral2;
integral1=0;
integral2=0;
for np=1:NP
    z_DB=del*np*1i-a*1i/2+a/2;
    L1=zeros(N,1);
    L2=zeros(M,1);
    for n=1:N
        L1(n,1)=n*(R1/z_DB)^(n+1);
    end
    for m=1:M
        L2(m,1)=Get_P_Diff(m,a,z_DB);
    end    
    integral1=integral1+imag(Xe2.'*S2.'*L2-1/R1*Xd2.'*S1.'*L1)*del;
    integral2=integral2+imag(Xe1.'*S2.'*L2-1/R1*Xd1.'*S1.'*L1)*del;
end
delta_cof(2,1)=integral1;
delta_cof(2,2)=integral2;
Stress_res(1,1)=S023_Divided_G0;
Stress_res(2,1)=S013_Divided_G0;
res=delta_cof\Stress_res;
delta1_Divided_a=res(1,1);
delta2_Divided_a=res(2,1);
change=zeros(2,2);
change(1,2)=1;
change(2,1)=1;
delta_cof=change*delta_cof;
end